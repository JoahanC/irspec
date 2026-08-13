"""Resolution-match and stitch JWST MRS sub-band cubes into one spectral cube.

The MRS sub-bands are delivered as separate ``*_s3d`` cubes, each with its own
spatial WCS, pixel scale, and -- critically -- its own *angular resolution*: the
PSF FWHM grows roughly linearly with wavelength, from ~0.3" in channel 1 (~5 um)
to ~1.0" in channel 4 (~28 um). A fixed spaxel therefore samples a different
physical region at every wavelength, which corrupts spectral maps and per-spaxel
SEDs.

:class:`CubeStitcher` brings every wavelength slice to a single common angular
resolution (the coarsest, i.e. the reddest band) and a single common spatial
grid, then concatenates the aligned slices along wavelength into one
``(n_wave, ny, nx)`` cube in which each spaxel is a continuous, resolution-matched
spectrum.

Pipeline (see :meth:`CubeStitcher.build`):

1. Choose the target angular resolution (coarsest band) and a common tangent-plane
   WCS covering the union footprint of all input cubes.
2. PSF-homogenise every slice: convolve with a Gaussian kernel that brings its
   PSF up to the target FWHM (skip the reddest band, which is already coarsest).
3. Reproject every (convolved) slice onto the common spatial grid.
4. Stitch: order all slices by wavelength, drop duplicated band-overlap channels,
   optionally apply per-band flux-matching ratios.

This is a first-pass draft. The transforms in steps 1-3 are functional; the
flux-matching in step 4 and rigorous error propagation are marked below as
``DRAFT`` and need refinement before science use (see caveats in each method).

Notes:
    The input cubes are assumed to be the raw ``s3d`` products in the *observed*
    frame (``redshift=0`` on the :class:`~irspec.datacube.Datacube`), so the PSF
    model is evaluated at observed wavelength.
"""

import re
import warnings
from dataclasses import dataclass

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.wcs.utils import proj_plane_pixel_scales

# JWST primary aperture, for the diffraction-limited fallback PSF model.
_JWST_DIAMETER = 6.5 * u.m
_FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))

# Pipeline calibration version (``CAL_VER``) at which the MIRI/MRS covariance-
# aware ERR scaling (Law et al. 2023) was incorporated -- operational Build 11.0,
# i.e. jwst 1.14.0. Cubes calibrated *before* this carry per-voxel ERR that both
# under-estimates the true per-voxel noise and ignores the large voxel-to-voxel
# covariance from cube resampling; ``err_rescale='auto'`` corrects only those.
_ERR_FIX_CALVER = (1, 14, 0)


def miri_mrs_fwhm(wave_um):
    """MIRI/MRS PSF FWHM in arcsec at ``wave_um`` (Law et al. 2023, ~linear)."""
    return 0.033 * np.asarray(wave_um, dtype=float) + 0.106


def diffraction_fwhm(wave_um):
    """Diffraction-limited FWHM (arcsec): ~1.03 lambda/D for a filled aperture."""
    wave = np.asarray(wave_um, dtype=float) * u.micron
    theta = (1.03 * wave / _JWST_DIAMETER).to(u.dimensionless_unscaled).value
    return np.degrees(theta) * 3600.0


@dataclass
class StitchedCube:
    """Result of :meth:`CubeStitcher.build`.

    Attributes:
        flux (ndarray): The stitched cube, shape ``(n_wave, ny, nx)``.
        error (ndarray): Propagated 1-sigma errors, same shape (see caveat in
            :meth:`CubeStitcher._reproject_slice`).
        dq (ndarray): Data-quality/coverage flags, same shape.
        waves (ndarray): Wavelengths in micron, length ``n_wave``, ascending.
        wcs (WCS): The 2-D celestial WCS of the common spatial grid.
        bandnames (ndarray): Source sub-band label per wavelength channel.
        stitch_log (list): Per-band record of the flux-matching correction
            applied (band name, whether it overlapped its bluer neighbour, and
            the correction summary). Empty when ``flux_match=False``.
        err_rescale_log (list): Per-band record of the input-ERR rescaling
            (band name, applied ``factor``, and a ``note``). Empty/factor 1.0
            when ``err_rescale`` is disabled or the cubes are post-fix.
    """
    flux: np.ndarray
    error: np.ndarray
    dq: np.ndarray
    waves: np.ndarray
    wcs: WCS
    bandnames: np.ndarray
    stitch_log: list = None
    err_rescale_log: list = None

    def write(self, path):
        """Write a multi-extension FITS cube (FLUX/ERR/DQ/WAVE) with the common
        3-D WCS, mirroring the CRETA ``*_cube.fits`` layout so downstream
        tooling can read it."""
        hdr = self.wcs.to_header()
        primary = fits.PrimaryHDU()
        # Record any input-ERR rescaling in the primary header for provenance.
        for e in (self.err_rescale_log or []):
            if e.get("factor", 1.0) != 1.0 or e.get("note"):
                primary.header.add_history(
                    f"ERR rescale {e['band']}: x{e.get('factor', 1.0):.3f}"
                    f" ({e.get('note', '')})")
        flux_hdu = fits.ImageHDU(self.flux, header=hdr, name="FLUX")
        hdus = [primary, flux_hdu,
                fits.ImageHDU(self.error, header=hdr, name="ERR"),
                fits.ImageHDU(self.dq.astype("float32"), header=hdr, name="DQ"),
                fits.ImageHDU(self.waves.astype("float64"), name="WAVE")]
        col = fits.Column(name="Band_name", format="20A", array=self.bandnames)
        hdus.append(fits.BinTableHDU.from_columns([col], name="Band_names"))
        fits.HDUList(hdus).writeto(path, overwrite=True)
        return path


class CubeStitcher:
    """Resolution-match a set of MRS sub-band cubes and stitch them into one
    spectral cube on a common spatial grid at a common angular resolution.
    """

    def __init__(self, cubes, psf_match=True, sampling_frac=0.5,
                 fwhm_model=miri_mrs_fwhm, pad_pix=2,
                 flux_match=True, match_deg=1, reference="snr",
                 per_spaxel=False, flux_method="chained",
                 joint_iters=30, joint_grid=2000,
                 hybrid_snr_lo=8.0, hybrid_snr_hi=40.0,
                 err_rescale=None):
        """Initialise the stitcher.

        Args:
            cubes (list[Datacube]): The sub-band cubes to combine, in any
                wavelength order. Each must expose ``science_data``,
                ``error_data``, ``dq_data``, ``wcs`` (2-D celestial), ``wvs``,
                ``im_shape`` and ``_sb_unit`` (the :class:`Datacube` interface).
            psf_match (bool): If ``True``, convolve each slice up to the common
                angular resolution before reprojecting. If ``False``, only
                resample onto the common grid (matches pixel sampling, not PSF).
            sampling_frac (float): Common pixel scale as a fraction of the target
                FWHM (0.5 = Nyquist). Smaller oversamples the coarse PSF.
            fwhm_model (callable): ``fwhm_arcsec = fwhm_model(wave_um)``; defaults
                to the MIRI/MRS relation. Use :func:`diffraction_fwhm` or a custom
                callable for other instruments.
            pad_pix (int): Extra pixels added around the union footprint.
            flux_match (bool): If ``True``, remove residual band-to-band flux
                jumps by anchoring every band to a reference band across their
                overlaps (see :meth:`_apply_flux_matching`).
            match_deg (int): Polynomial degree of the per-band flux-matching
                correction across the overlap (0 = scalar, 1 = linear, ...).
            reference (int | str): Which band anchors the flux scale:
                ``'snr'`` (default, highest median SNR -- the most reliable
                anchor), ``'blue'`` / ``'red'`` (bluest/reddest band), a band
                name (e.g. ``'ch2_SHORT'``), or an integer index. Bands are
                corrected outward from the reference in both directions.
            per_spaxel (bool): Experimental (default ``False``). Modulate each
                band's flux-matching correction by a *gentle* per-spaxel scale
                (SNR-gated, spatially smoothed, bounded to +/-15%) for residual
                spatially-varying offsets. Off by default because the overlap
                regions are the noisy band roll-off edges, where per-spaxel
                ratios are unreliable; the dominant band-edge jumps are handled
                instead by the midpoint-crossover trim (see
                :meth:`_band_keep_masks`).
            flux_method (str): How the band scales are estimated.
                ``'chained'`` (default) walks outward from the reference, tying
                each band to its already-fixed neighbour over their overlap
                (see :meth:`_apply_flux_matching`). ``'joint'`` instead solves a
                rank-1 matrix-completion: one shared spectral template and one
                scale per band (per spaxel when ``per_spaxel``), so every band is
                tied to the *global* SED rather than only its neighbour. On MRS's
                path-shaped overlap graph the two agree for clean data, but the
                joint fit does not accumulate error along the chain, so it is more
                robust for faint/noisy spaxels (see :meth:`_apply_joint_rank1`).
                ``'hybrid'`` blends the two: it follows the low-variance local
                chain where the overlap SNR is high and falls back to the joint
                scale where the overlap is faint, giving the best of both on faint
                data (see :meth:`_apply_hybrid`).
            joint_iters (int): Alternating-update iterations for the ``'joint'``
                rank-1 solve.
            joint_grid (int): Number of points on the common log-wavelength grid
                the ``'joint'`` template is built on.
            hybrid_snr_lo (float): Overlap SNR at/below which the ``'hybrid'``
                method fully trusts the joint fallback (gate weight 0).
            hybrid_snr_hi (float): Overlap SNR at/above which the ``'hybrid'``
                method fully trusts the local chained ratio (gate weight 1).
            err_rescale (None | float | str): Optionally inflate each input
                cube's ``ERR`` before stitching, to compensate for the known
                MIRI/MRS + NIRSpec IFU pipeline ERR under-estimation. ``None``
                (default) leaves errors untouched. A float multiplies every
                band's ERR by that fixed factor. ``'auto'`` measures each band's
                empirical-to-pipeline noise ratio ``R1`` (robust spectral
                2nd-difference noise / median ERR) and applies ``max(R1, 1)``,
                clipped to ``[1, 5]`` -- but **only for cubes calibrated before
                the Build 11.0 covariance fix** (``CAL_VER < 1.14.0``); post-fix
                or unrecognised cubes are left unchanged, so the same call is a
                no-op on reprocessed data. The rescaled ERR then feeds the anchor
                choice, the hybrid overlap-SNR gate, and the propagated errors.
                See :meth:`_err_rescale_factor`.
        """
        if len(cubes) == 0:
            raise ValueError("No cubes provided to stitch.")
        # Order by ascending starting wavelength (fix #2): the sequential
        # flux-matching and overlap-trimming both rely on wavelength order.
        self.cubes = sorted(cubes, key=lambda c: float(np.min(c.wvs.value)))
        self.psf_match = psf_match
        self.sampling_frac = sampling_frac
        self.fwhm_model = fwhm_model
        self.pad_pix = pad_pix
        self.flux_match = flux_match
        self.match_deg = match_deg
        self.reference = reference
        self.per_spaxel = per_spaxel
        self.flux_method = flux_method
        self.joint_iters = joint_iters
        self.joint_grid = joint_grid
        self.hybrid_snr_lo = hybrid_snr_lo
        self.hybrid_snr_hi = hybrid_snr_hi
        self.err_rescale = err_rescale
        self.stitch_log = []
        self.err_rescale_log = []

    # -- resolution + grid ---------------------------------------------------

    def _sigma_arcsec(self, wave_um):
        """PSF sigma (arcsec) at ``wave_um`` from the FWHM model."""
        return self.fwhm_model(wave_um) / _FWHM_PER_SIGMA

    def _target_resolution(self):
        """The coarsest angular resolution across all bands, as an arcsec FWHM
        (evaluated at each band's reddest wavelength, which is its coarsest)."""
        return max(float(self.fwhm_model(np.max(c.wvs.value))) for c in self.cubes)

    def _common_wcs(self):
        """Build the common tangent-plane WCS + output shape covering the union
        footprint of all cubes at ``sampling_frac * target_FWHM`` per pixel.

        Returns:
            tuple: ``(wcs, ny, nx)`` for the common grid.
        """
        target_fwhm = self._target_resolution()
        pix_arcsec = self.sampling_frac * target_fwhm
        pix_deg = pix_arcsec / 3600.0

        # Collect the sky corners of every cube.
        ras, decs = [], []
        for cube in self.cubes:
            ny, nx = cube.im_shape
            xs = [0, nx - 1, 0, nx - 1]
            ys = [0, 0, ny - 1, ny - 1]
            sky = cube.wcs.pixel_to_world(xs, ys)
            ras.append(np.asarray(sky.ra.deg))
            decs.append(np.asarray(sky.dec.deg))
        ra = np.concatenate(ras)
        dec = np.concatenate(decs)

        center = SkyCoord(ra.mean() * u.deg, dec.mean() * u.deg)
        cos_dec = np.cos(np.radians(center.dec.deg))
        # Extent on the sky (arcsec), RA scaled by cos(dec).
        span_ra = (ra.max() - ra.min()) * cos_dec
        span_dec = dec.max() - dec.min()
        nx = int(np.ceil(span_ra / pix_deg)) + 2 * self.pad_pix
        ny = int(np.ceil(span_dec / pix_deg)) + 2 * self.pad_pix

        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crval = [center.ra.deg, center.dec.deg]
        wcs.wcs.crpix = [(nx + 1) / 2.0, (ny + 1) / 2.0]  # 1-indexed centre
        wcs.wcs.cdelt = [-pix_deg, pix_deg]  # RA increases to the left
        return wcs, ny, nx

    # -- per-slice transforms ------------------------------------------------

    def _homogenize_slice(self, slice2d, pix_scale_arcsec, sigma_band, sigma_tgt,
                          variance=False):
        """Convolve a slice up to the target resolution.

        Uses a Gaussian kernel of width ``sqrt(sigma_tgt**2 - sigma_band**2)``.
        If the band is already at/above the target resolution the slice is
        returned unchanged (you can only broaden, never sharpen).

        Args:
            slice2d (ndarray): 2-D image (may contain NaNs at FOV edges).
            pix_scale_arcsec (float): Pixel scale of this slice, arcsec/pix.
            sigma_band (float): This slice's PSF sigma, arcsec.
            sigma_tgt (float): Target PSF sigma, arcsec.
            variance (bool): If ``True``, treat ``slice2d`` as a variance map
                and convolve with the *squared* (un-normalised) kernel, so the
                result is ``Var(sum_i w_i x_i) = sum_i w_i**2 Var(x_i)`` for
                independent pixels. Pass ``err**2`` and take ``sqrt`` afterwards.
                (Still neglects pre-existing inter-pixel covariance.)

        Returns:
            ndarray: The (possibly) convolved slice.
        """
        if sigma_tgt <= sigma_band:
            return slice2d
        kernel_sigma_arcsec = np.sqrt(sigma_tgt**2 - sigma_band**2)
        kernel_sigma_pix = kernel_sigma_arcsec / pix_scale_arcsec
        if kernel_sigma_pix < 1e-3:
            return slice2d
        kernel = Gaussian2DKernel(kernel_sigma_pix)
        if variance:
            # Propagate variance: convolve with the squared weights (sum != 1).
            from astropy.convolution import CustomKernel
            k = CustomKernel(kernel.array ** 2)
            return convolve_fft(slice2d, k, nan_treatment="interpolate",
                                preserve_nan=True, normalize_kernel=False,
                                allow_huge=True)
        # NaN-aware, flux-preserving (kernel is normalised).
        return convolve_fft(slice2d, kernel, nan_treatment="interpolate",
                            preserve_nan=True, normalize_kernel=True,
                            allow_huge=True)

    def _reproject_slice(self, slice2d, wcs_in, wcs_out, shape_out, flux_per_pixel):
        """Resample a slice onto the common grid.

        Args:
            slice2d (ndarray): 2-D image to resample.
            wcs_in (WCS): Input celestial WCS.
            wcs_out (WCS): Target celestial WCS.
            shape_out (tuple): ``(ny, nx)`` of the target grid.
            flux_per_pixel (bool): ``True`` if the data are per-pixel flux
                (e.g. Jy/pix), which must be flux-conserved; ``False`` for
                surface brightness (e.g. MJy/sr), which interpolates.

        Returns:
            ndarray: The resampled slice on the common grid.

        DRAFT caveat: interpolation/convolution correlate neighbouring pixels,
        so a variance array pushed through the same path *underestimates* the
        true (now-covariant) uncertainty. A rigorous version must track the
        covariance or inflate the errors accordingly.
        """
        if flux_per_pixel:
            from reproject import reproject_exact
            out, _ = reproject_exact((slice2d, wcs_in), wcs_out,
                                     shape_out=shape_out, parallel=False)
        else:
            from reproject import reproject_interp
            out, _ = reproject_interp((slice2d, wcs_in), wcs_out,
                                      shape_out=shape_out, order="bilinear")
        return out

    # -- stitching -----------------------------------------------------------

    @staticmethod
    def _integrated_spectrum(band):
        """Spatially-averaged spectrum (nanmean over all spaxels) of a band."""
        with np.errstate(invalid="ignore"):
            spec = np.nanmean(band["flux"], axis=(1, 2))
        return band["waves"], spec

    @staticmethod
    def _coverage_mask(band, lo, hi):
        """Spaxels a band covers over ``[lo, hi]`` (finite in >50% of the
        overlap channels)."""
        m = (band["waves"] >= lo) & (band["waves"] <= hi)
        if not m.any():
            m = np.ones(len(band["waves"]), bool)
        with np.errstate(invalid="ignore"):
            frac = np.isfinite(band["flux"][m]).mean(axis=0)
        return frac > 0.5

    @staticmethod
    def _masked_spectrum(band, mask):
        """Mean spectrum over a fixed set of spaxels (unbiased by footprint)."""
        with np.errstate(invalid="ignore"):
            return np.nanmean(band["flux"][:, mask], axis=1)

    @staticmethod
    def _masked_error(band, mask):
        """Mean per-channel error over a fixed set of spaxels."""
        with np.errstate(invalid="ignore"):
            return np.sqrt(np.nanmean(band["err"][:, mask] ** 2, axis=1))

    @staticmethod
    def _band_snr(band):
        """Median per-pixel signal-to-noise of a band, used to pick the anchor."""
        with np.errstate(invalid="ignore", divide="ignore"):
            snr = np.where(band["err"] > 0, band["flux"] / band["err"], np.nan)
        return float(np.nanmedian(snr))

    # -- input ERR rescaling (pipeline ERR under-estimation) -----------------

    @staticmethod
    def _calver(cube):
        """The cube's pipeline ``CAL_VER`` as an ``(major, minor, patch)`` tuple,
        or ``None`` if absent/unparseable."""
        hdr = getattr(cube, "general_header", None)
        val = hdr.get("CAL_VER") if hdr else None
        if not val:
            return None
        m = re.match(r"\s*(\d+)\.(\d+)\.(\d+)", str(val))
        return tuple(int(g) for g in m.groups()) if m else None

    @staticmethod
    def _empirical_noise(spec):
        """Robust 1-sigma noise of a 1-D spectrum from its spectral 2nd
        difference (``f[i-1] - 2 f[i] + f[i+1]``), which suppresses smooth
        signal so the estimate reflects pixel noise, not source structure.
        Scaled by ``1/sqrt(6)`` for the 2nd-difference variance inflation and by
        1.482 (MAD -> sigma)."""
        d2 = spec[2:] - 2.0 * spec[1:-1] + spec[:-2]
        d2 = d2[np.isfinite(d2)]
        if d2.size < 8:
            return np.nan
        return 1.482 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6.0)

    def _measure_err_ratio(self, cube):
        """Median over the brighter half of the FOV of ``empirical_noise / ERR``
        for a cube -- the factor by which the pipeline ERR under-estimates the
        true per-voxel noise (``R1`` in the covariance test). Measured on the
        native (pre-resample) product so the resampling covariance does not
        contaminate it; a scalar factor commutes with the later PSF/reproject."""
        sci = np.asarray(cube.science_data, dtype=float)
        err = np.asarray(cube.error_data, dtype=float)
        nw = sci.shape[0]
        lo, hi = int(0.20 * nw), int(0.80 * nw)      # trim noisy band roll-offs
        if hi - lo < 20:
            lo, hi = 0, nw
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            med_img = np.nanmedian(sci[lo:hi], axis=0)
            emed = np.nanmedian(err[lo:hi], axis=0)
            finite = np.isfinite(med_img) & (emed > 0)
            if not finite.any():
                return np.nan
            thr = np.nanmedian(med_img[finite])
        ys, xs = np.where(finite & (med_img > thr))
        if ys.size < 4:                              # too few bright spaxels
            ys, xs = np.where(finite)
        ratios = []
        for y, x in zip(ys, xs):
            sn = self._empirical_noise(sci[lo:hi, y, x])
            e = emed[y, x]
            if np.isfinite(sn) and e > 0:
                ratios.append(sn / e)
        return float(np.nanmedian(ratios)) if ratios else np.nan

    def _err_rescale_factor(self, cube):
        """Resolve ``self.err_rescale`` to a ``(factor, note)`` for one cube.

        ``None`` -> ``(1.0, None)`` (disabled). A number -> that fixed factor,
        applied unconditionally. ``'auto'`` -> build-gated empirical rescale:
        only cubes calibrated before :data:`_ERR_FIX_CALVER` (the Build 11.0
        covariance fix) are rescaled by ``clip(max(R1, 1), 1, 5)``; post-fix or
        unrecognised cubes get ``1.0`` with an explanatory note.
        """
        mode = self.err_rescale
        band = self._bandname(cube)
        if mode is None or mode is False:
            return 1.0, None
        if isinstance(mode, bool):                   # True is not a valid factor
            raise ValueError("err_rescale must be None, a float, or 'auto'.")
        if isinstance(mode, (int, float)):
            return float(mode), f"forced x{float(mode):.2f}"
        if mode == "auto":
            calver = self._calver(cube)
            vstr = ".".join(map(str, calver)) if calver else "unknown"
            if calver is None or calver >= _ERR_FIX_CALVER:
                return 1.0, (f"auto: CAL_VER {vstr} is post-fix/unknown "
                             f"-> no rescale")
            r = self._measure_err_ratio(cube)
            if not np.isfinite(r):
                return 1.0, (f"auto: pre-fix (CAL_VER {vstr}) but R1 "
                             f"unmeasurable -> no rescale")
            factor = float(np.clip(max(r, 1.0), 1.0, 5.0))
            return factor, (f"auto: pre-fix (CAL_VER {vstr}), R1={r:.2f} "
                            f"-> ERR x{factor:.2f}")
        raise ValueError(f"Unknown err_rescale={mode!r}; use None, a float, "
                         f"or 'auto'.")

    def _reference_index(self, band_cubes):
        """Resolve ``self.reference`` to a band index: an int, a band name, or
        ``'blue'`` / ``'red'`` / ``'snr'`` (highest median SNR, the default and
        most reliable anchor)."""
        ref = self.reference
        if isinstance(ref, (int, np.integer)):
            return int(ref) % len(band_cubes)
        names = [b["name"] for b in band_cubes]
        if ref in names:
            return names.index(ref)
        if ref == "blue":
            return 0
        if ref == "red":
            return len(band_cubes) - 1
        if ref == "snr":
            return int(np.argmax([self._band_snr(b) for b in band_cubes]))
        raise ValueError(f"Unknown reference {ref!r}; use an int, a band name, "
                         f"or 'blue'/'red'/'snr'.")

    def _match_pair(self, anchor, target):
        """Scale ``target`` (in place) onto ``anchor``'s already-fixed flux
        scale via a flux-matching correction measured on their overlap.
        Returns a log entry.

        * **Linear-over-overlap (#3):** fit the per-wavelength flux ratio of the
          common-footprint spectra with a degree-``match_deg`` polynomial
          (removes a tilt, not just an offset), held constant outside the
          overlap window. This sets the global wavelength shape ``C(lambda)``.
        * **Per-spaxel modulation (refinement):** when ``per_spaxel`` is set,
          ``C(lambda)`` is multiplied by a per-spaxel scale ``k(y,x)`` -- the
          ratio of the two bands' overlap-integrated flux at that spaxel,
          normalised so its footprint median is 1 -- so *every* spaxel is made
          continuous, not just the spatial average. The global normalisation
          (and hence the anchoring) is unchanged.
        * **Error propagation (refinement):** the correction's own fractional
          uncertainty (the scatter of the overlap ratio about the fit) is folded
          into the matched errors.
        """
        lo = max(float(anchor["waves"].min()), float(target["waves"].min()))
        hi = min(float(anchor["waves"].max()), float(target["waves"].max()))
        common = (self._coverage_mask(anchor, lo, hi)
                  & self._coverage_mask(target, lo, hi)) if hi > lo else None
        if common is None or int(common.sum()) < 3:
            return {"band": target["name"], "overlap": False,
                    "correction": "no usable overlap -> left on native scale"}

        a_w, a_s = anchor["waves"], self._masked_spectrum(anchor, common)
        t_w, t_s = target["waves"], self._masked_spectrum(target, common)
        af = np.isfinite(a_s)
        mo = (t_w >= lo) & (t_w <= hi)
        a_interp = np.interp(t_w[mo], a_w[af], a_s[af])
        t_over = t_s[mo]
        ratio = np.divide(a_interp, t_over, out=np.full_like(t_over, np.nan),
                          where=t_over > 0)
        good = (np.isfinite(ratio) & (a_interp > 0)
                & (ratio > 0.1) & (ratio < 10.0))  # reject line/edge spikes
        if good.sum() < 2:
            coeffs, rel_sigma = np.array([1.0]), 0.0  # no reliable correction
        else:
            deg = self.match_deg if good.sum() > self.match_deg + 1 else 0
            # SNR-weight the fit so the noisy band roll-off channels in the
            # overlap do not drive the correction.
            t_e_over = self._masked_error(target, common)[mo]
            weights = np.where(t_e_over[good] > 0, t_over[good] / t_e_over[good], 1.0)
            coeffs = np.polyfit(t_w[mo][good], ratio[good], deg, w=weights)
            resid = ratio[good] - np.polyval(coeffs, t_w[mo][good])
            rel_sigma = float(np.clip(np.std(resid) / np.mean(ratio[good]), 0.0, 0.5))

        # Global wavelength shape C(lambda), held constant outside the overlap.
        corr = np.polyval(coeffs, t_w)
        corr[t_w < lo] = float(np.polyval(coeffs, lo))
        corr[t_w > hi] = float(np.polyval(coeffs, hi))
        corr = np.clip(corr, 0.2, 5.0)  # guard against pathological fits
        corr3d = corr[:, None, None]

        # Per-spaxel modulation k(y,x): ratio of the bands' overlap-integrated
        # flux per spaxel, normalised to a footprint median of 1.
        k_note = ""
        if self.per_spaxel:
            from scipy.ndimage import gaussian_filter
            a_mo = (a_w >= lo) & (a_w <= hi)
            with np.errstate(invalid="ignore", divide="ignore"):
                a_map = np.nanmean(anchor["flux"][a_mo], axis=0)
                t_map = np.nanmean(target["flux"][mo], axis=0)
                t_emap = np.sqrt(np.nanmean(target["err"][mo] ** 2, axis=0))
                s = np.divide(a_map, t_map, out=np.full_like(t_map, np.nan),
                              where=t_map > 0)
                tsnr = np.where(t_emap > 0, t_map / t_emap, 0.0)
            s_ref = np.nanmedian(s[common])
            k = np.where(np.isfinite(s) & (s > 0) & (s_ref > 0), s / s_ref, 1.0)
            # SNR-gate: shrink toward 1 where the overlap is faint/noisy.
            w = np.clip((tsnr - 5.0) / 15.0, 0.0, 1.0)
            k = 1.0 + (k - 1.0) * w
            k = gaussian_filter(np.nan_to_num(k, nan=1.0), sigma=0.8)  # de-noise
            k = np.clip(k, 0.85, 1.18)  # a gentle refinement, never a big swing
            corr3d = corr3d * k[None, :, :]
            k_note = f", per-spaxel k[{k[common].min():.2f},{k[common].max():.2f}]"

        f0 = target["flux"]
        target["flux"] = f0 * corr3d
        target["err"] = np.abs(corr3d) * np.sqrt(target["err"]**2
                                                 + (rel_sigma * f0)**2)
        return {"band": target["name"], "overlap": True,
                "correction": f"poly deg{len(coeffs) - 1} over [{lo:.2f},{hi:.2f}] um "
                              f"({int(common.sum())} common spaxels), factor "
                              f"{float(np.polyval(coeffs, lo)):.3f}"
                              f"->{float(np.polyval(coeffs, hi)):.3f}"
                              f", err+{rel_sigma * 100:.0f}%{k_note}"}

    def _apply_flux_matching(self, band_cubes):
        """Anchor every band onto a chosen reference band's flux scale, in place.

        Improvements over CRETA's single scalar median ratio:

        * **Configurable reference (#4).** The anchor band is chosen by
          ``reference`` (default ``'snr'``, the highest-SNR band) instead of
          always the bluest, so a possibly-anomalous edge band's calibration is
          not propagated to the whole SED. Bands are corrected *outward* from the
          reference in both directions (redward and blueward).
        * **Linear-over-overlap correction (#3)** on the **common spatial
          footprint** -- see :meth:`_match_pair`.
        * **No de-anchoring on missing overlap (#1).** A pair with no overlap
          leaves that band native and the chain continues; it never rescales the
          already-matched bands (unlike CRETA's ``1.0`` sentinel).

        The single per-band correction is applied to every spaxel (a per-spaxel
        correction is a future refinement). Mutates ``band_cubes``; returns a
        per-band log in wavelength order.
        """
        r = self._reference_index(band_cubes)
        log = [None] * len(band_cubes)
        log[r] = {"band": band_cubes[r]["name"], "overlap": None,
                  "correction": f"reference / anchor ({self.reference})"}
        # Redward of the reference: each band tied to its bluer neighbour.
        for k in range(r + 1, len(band_cubes)):
            log[k] = self._match_pair(band_cubes[k - 1], band_cubes[k])
        # Blueward of the reference: each band tied to its redder neighbour.
        for k in range(r - 1, -1, -1):
            log[k] = self._match_pair(band_cubes[k + 1], band_cubes[k])
        return log

    # -- joint rank-1 flux matching (Haute-Couture-style) --------------------

    def _common_loggrid(self, band_cubes):
        """A common log-spaced wavelength grid spanning all bands, on which the
        shared rank-1 template is built."""
        lo = min(float(b["waves"][0]) for b in band_cubes)
        hi = max(float(b["waves"][-1]) for b in band_cubes)
        return np.geomspace(lo, hi, int(self.joint_grid))

    @staticmethod
    def _resample_band(waves, arr, grid):
        """Linearly resample ``arr`` (axis-0 = wavelength ``waves``, monotonic)
        onto ``grid`` along axis 0, returning NaN outside the band's support.
        Vectorised over any trailing (spatial) axes."""
        j = np.clip(np.searchsorted(waves, grid) - 1, 0, len(waves) - 2)
        with np.errstate(invalid="ignore"):
            w = (grid - waves[j]) / (waves[j + 1] - waves[j])
        w = w.reshape((-1,) + (1,) * (arr.ndim - 1))
        out = (1.0 - w) * arr[j] + w * arr[j + 1]
        out[(grid < waves[0]) | (grid > waves[-1])] = np.nan
        return out

    @staticmethod
    def _rank1_scales(Fg, Wg, ref, iters=30):
        """Rank-1 matrix completion for per-band flux scales.

        Alternately (i) build one shared spectral template ``T`` as the
        weight-agnostic median across bands of the currently-scaled fluxes, and
        (ii) update each band's scale to the SNR-weighted geometric mean of
        ``T / flux`` over that band's support. Gauge-fixed so the reference
        band's scale is 1. Because every band is matched to the *global* template
        rather than only its neighbour, the estimate does not accumulate error
        along the band chain.

        Args:
            Fg (ndarray): Fluxes on the common grid, shape ``(n_band, n_grid,
                *spatial)``; NaN where a band does not cover a grid point.
            Wg (ndarray): Non-negative weights (per-point SNR), same shape.
            ref (int): Reference band index (its scale is pinned to 1).
            iters (int): Alternating-update iterations.

        Returns:
            ndarray: Scales, shape ``(n_band, *spatial)``.
        """
        valid = np.isfinite(Fg) & (Fg > 0)
        Fpos = np.where(valid, Fg, np.nan)
        w = np.where(valid, np.nan_to_num(np.clip(Wg, 0.0, None), nan=0.0), 0.0)
        s = np.ones((Fg.shape[0],) + Fg.shape[2:])
        for _ in range(iters):
            with np.errstate(invalid="ignore"), warnings.catch_warnings():
                # grid points outside every band's coverage give All-NaN medians.
                warnings.simplefilter("ignore", category=RuntimeWarning)
                T = np.nanmedian(Fpos * s[:, None], axis=0)             # template
                lr = np.log(np.clip(T[None] / Fpos, 0.1, 10.0))
            num = np.nansum(w * lr, axis=1)
            den = np.nansum(w, axis=1)
            s = np.exp(np.divide(num, den, out=np.zeros_like(num), where=den > 0))
            sref = s[ref]
            s = np.divide(s, sref[None], out=np.ones_like(s), where=sref[None] > 0)
            s = np.where(np.isfinite(s) & (s > 0), s, 1.0)
        return s

    def _joint_scales(self, band_cubes, ref):
        """Rank-1 joint scales for every band, plus a descriptive note.

        Builds a common log-wavelength grid, resamples every band's flux and SNR
        onto it, and solves :meth:`_rank1_scales`. Always returns the global scale
        (from the footprint-averaged spectra); when ``per_spaxel`` is set it also
        runs the per-spaxel solve and *shrinks it toward the global scale* on
        faint spaxels (gated by each spaxel's median SNR), so noisy spaxels fall
        back to the robust global correction instead of chasing noise.

        Returns:
            tuple: ``(scales, note)`` where ``scales`` is a list -- one float per
            band (global) or one ``(ny, nx)`` array per band (per-spaxel).
        """
        grid = self._common_loggrid(band_cubes)
        Fg = np.stack([self._resample_band(b["waves"], b["flux"], grid)
                       for b in band_cubes])
        Eg = np.stack([self._resample_band(b["waves"], b["err"], grid)
                       for b in band_cubes])
        with np.errstate(invalid="ignore", divide="ignore"):
            Wg = np.where(Eg > 0, Fg / Eg, 0.0)

        if Fg.ndim > 2:
            spatial_axes = tuple(range(2, Fg.ndim))
            with np.errstate(invalid="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                Fg_mean = np.nanmean(Fg, axis=spatial_axes)
                Wg_mean = np.nanmean(Wg, axis=spatial_axes)
        else:
            Fg_mean, Wg_mean = Fg, Wg
        s_glob = self._rank1_scales(Fg_mean, Wg_mean, ref, self.joint_iters)

        if self.per_spaxel and Fg.ndim > 2:
            s = self._rank1_scales(Fg, Wg, ref, self.joint_iters)      # (Nb,ny,nx)
            g = s_glob.reshape((len(band_cubes),) + (1,) * len(Fg.shape[2:]))
            with np.errstate(invalid="ignore"):
                snr_spx = np.nanmedian(np.where(np.isfinite(Wg), Wg, np.nan),
                                       axis=(0, 1))
            gate = np.clip((np.nan_to_num(snr_spx) - 3.0) / 12.0, 0.0, 1.0)
            s = g * (np.divide(s, g, out=np.ones_like(s), where=g > 0)
                     ** gate[None, ...])
            return [s[i] for i in range(len(band_cubes))], \
                "per-spaxel rank-1 (shrunk to global on faint spaxels)"
        return [float(s_glob[i]) for i in range(len(band_cubes))], "global rank-1"

    def _apply_scales(self, band_cubes, ref, scales, note, tag):
        """Apply one multiplicative scale per band (float or ``(ny,nx)`` map) to
        flux and error in place, and build the per-band log."""
        log = [None] * len(band_cubes)
        for i, b in enumerate(band_cubes):
            corr = scales[i]
            corr3d = corr[None, :, :] if np.ndim(corr) else corr
            b["flux"] = b["flux"] * corr3d
            b["err"] = b["err"] * np.abs(corr3d)
            if i == ref:
                log[i] = {"band": b["name"], "overlap": None,
                          "correction": f"reference / anchor ({self.reference}), {note}"}
            else:
                gmean = float(np.nanmean(corr)) if np.ndim(corr) else corr
                log[i] = {"band": b["name"], "overlap": True,
                          "correction": f"{tag} scale {gmean:.3f} ({note})"}
        return log

    def _apply_joint_rank1(self, band_cubes):
        """Flux-match all bands with the rank-1 joint solve (in place)."""
        ref = self._reference_index(band_cubes)
        scales, note = self._joint_scales(band_cubes, ref)
        return self._apply_scales(band_cubes, ref, scales, note, "joint rank-1")

    # -- hybrid: joint tail robustness + local override ----------------------

    def _hybrid_gate(self, snr):
        """Local-override weight in ``[0, 1]`` from an overlap SNR: 1 where the
        overlap is clean (trust the low-variance local ratio), 0 where it is
        faint (fall back to the robust joint scale)."""
        lo, hi = self.hybrid_snr_lo, self.hybrid_snr_hi
        return np.clip((snr - lo) / (hi - lo), 0.0, 1.0)

    def _overlap_scalar_ratio(self, anchor, target):
        """Robust scalar flux ratio (anchor/target) over their common-footprint
        overlap and the target's overlap SNR, or ``None`` if there is no usable
        overlap. Mirrors the degree-0 :meth:`_match_pair` estimate."""
        lo = max(float(anchor["waves"].min()), float(target["waves"].min()))
        hi = min(float(anchor["waves"].max()), float(target["waves"].max()))
        common = (self._coverage_mask(anchor, lo, hi)
                  & self._coverage_mask(target, lo, hi)) if hi > lo else None
        if common is None or int(common.sum()) < 3:
            return None
        a_w, a_s = anchor["waves"], self._masked_spectrum(anchor, common)
        t_w, t_s = target["waves"], self._masked_spectrum(target, common)
        af = np.isfinite(a_s)
        mo = (t_w >= lo) & (t_w <= hi)
        a_interp = np.interp(t_w[mo], a_w[af], a_s[af])
        t_over, t_e = t_s[mo], self._masked_error(target, common)[mo]
        ratio = np.divide(a_interp, t_over, out=np.full_like(t_over, np.nan),
                          where=t_over > 0)
        good = (np.isfinite(ratio) & (a_interp > 0)
                & (ratio > 0.1) & (ratio < 10.0) & (t_e > 0))
        if good.sum() < 1:
            return None
        w = t_over[good] / t_e[good]
        r = float(np.clip(np.sum(w * ratio[good]) / np.sum(w), 0.2, 5.0))
        return r, float(np.median(w))

    @staticmethod
    def _overlap_ratio_maps(anchor, target):
        """Per-spaxel (anchor/target) flux-ratio map over the overlap and the
        target's per-spaxel overlap-SNR map. Falls back to (1, 0) where a spaxel
        has no usable overlap, so the gate then selects the joint scale there."""
        lo = max(float(anchor["waves"].min()), float(target["waves"].min()))
        hi = min(float(anchor["waves"].max()), float(target["waves"].max()))
        shape = target["flux"].shape[1:]
        if hi <= lo:
            return np.ones(shape), np.zeros(shape)
        a_mo = (anchor["waves"] >= lo) & (anchor["waves"] <= hi)
        t_mo = (target["waves"] >= lo) & (target["waves"] <= hi)
        with np.errstate(invalid="ignore", divide="ignore"), \
                warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            a_map = np.nanmean(anchor["flux"][a_mo], axis=0)
            t_map = np.nanmean(target["flux"][t_mo], axis=0)
            t_emap = np.sqrt(np.nanmean(target["err"][t_mo] ** 2, axis=0))
            r = np.divide(a_map, t_map, out=np.full_like(t_map, np.nan),
                          where=t_map > 0)
            snr = np.where(t_emap > 0, t_map / t_emap, 0.0)
        r = np.clip(np.nan_to_num(r, nan=1.0), 0.2, 5.0)
        return r, np.nan_to_num(snr, nan=0.0)

    def _apply_hybrid(self, band_cubes):
        """Flux-match with the hybrid of joint rank-1 and local chaining (in
        place). Walking outward from the anchor, each band's scale is a geometric
        blend ``s_joint**(1-w) * (s_hybrid[nbr] * r)**w`` where ``r`` is the
        robust local overlap ratio to the already-fixed neighbour and ``w`` is the
        overlap-SNR gate. Clean overlaps follow the low-variance local chain;
        faint overlaps fall back to joint's absolute scale, which both avoids the
        chain's division-by-near-zero blow-up and stops it propagating outward.
        """
        ref = self._reference_index(band_cubes)
        s_joint, jnote = self._joint_scales(band_cubes, ref)
        per_spaxel = np.ndim(s_joint[ref]) > 0
        s_hyb = [None] * len(band_cubes)
        s_hyb[ref] = np.ones_like(s_joint[ref]) if per_spaxel else 1.0

        def blend(b, nbr):
            if per_spaxel:
                r, snr = self._overlap_ratio_maps(band_cubes[nbr], band_cubes[b])
            else:
                res = self._overlap_scalar_ratio(band_cubes[nbr], band_cubes[b])
                if res is None:
                    s_hyb[b] = s_joint[b]      # no overlap -> joint absolute
                    return
                r, snr = res
            w = self._hybrid_gate(snr)
            loc = np.clip(s_hyb[nbr] * r, 1e-6, None)
            s_hyb[b] = s_joint[b] ** (1.0 - w) * loc ** w

        for b in range(ref + 1, len(band_cubes)):
            blend(b, b - 1)
        for b in range(ref - 1, -1, -1):
            blend(b, b + 1)
        note = f"hybrid: {jnote} + local override (gate {self.hybrid_snr_lo:g}-" \
               f"{self.hybrid_snr_hi:g})"
        return self._apply_scales(band_cubes, ref, s_hyb, note, "hybrid")

    @staticmethod
    def _band_keep_masks(band_cubes):
        """Per-band channel masks that join adjacent bands at the **midpoint of
        their overlap** rather than keep-first.

        MRS sub-bands roll off (low throughput, rising noise) at their extreme
        edges. Keep-first retains the bluer band's declining red roll-off right
        up against the redder band's blue edge, which shows up as a flux jump at
        the join. Switching bands at the overlap midpoint drops *both* roll-off
        edges, landing the boundary where both bands are cleanest. Returns one
        boolean mask per band; the concatenation stays monotonic in wavelength.
        """
        masks = [np.ones(len(b["waves"]), bool) for b in band_cubes]
        for i in range(len(band_cubes) - 1):
            a, b = band_cubes[i], band_cubes[i + 1]
            lo = max(float(a["waves"].min()), float(b["waves"].min()))
            hi = min(float(a["waves"].max()), float(b["waves"].max()))
            if hi <= lo:
                continue  # no overlap -> keep both fully (a gap remains)
            mid = 0.5 * (lo + hi)
            masks[i] &= a["waves"] <= mid       # drop bluer band's red roll-off
            masks[i + 1] &= b["waves"] > mid     # drop redder band's blue roll-off
        return masks

    def build(self):
        """Run the full resolution-match + stitch pipeline.

        Returns:
            StitchedCube: The resolution-matched, spatially-aligned, stitched
            cube on the common grid.
        """
        wcs_out, ny, nx = self._common_wcs()
        sigma_tgt = self._target_resolution() / _FWHM_PER_SIGMA

        # 1. Resolution-match + reproject each band onto the common grid, kept
        #    separate so the overlap regions survive for flux-matching. Bands are
        #    in ascending-wavelength order (set in __init__, fix #2).
        self.err_rescale_log = []
        band_cubes = []
        for cube in self.cubes:
            name = self._bandname(cube)
            pix_scale = float(np.mean(proj_plane_pixel_scales(cube.wcs)) * 3600.0)
            flux_per_pixel = "pix" in str(getattr(cube, "_sb_unit", "")).lower()
            # Optionally inflate the (under-estimated) pipeline ERR before it
            # feeds the resolution match, flux match, and SNR-driven decisions.
            efac, enote = self._err_rescale_factor(cube)
            self.err_rescale_log.append(
                {"band": name, "factor": efac, "note": enote})
            fk, ek, wk = [], [], []
            for k in range(cube.ax3_len):
                wl = float(cube.wvs.value[k])
                fsl = cube.science_data[k]
                var = (cube.error_data[k] * efac) ** 2
                if self.psf_match:
                    sig_b = self._sigma_arcsec(wl)
                    fsl = self._homogenize_slice(fsl, pix_scale, sig_b, sigma_tgt)
                    var = self._homogenize_slice(var, pix_scale, sig_b, sigma_tgt,
                                                 variance=True)
                fk.append(self._reproject_slice(
                    fsl, cube.wcs, wcs_out, (ny, nx), flux_per_pixel))
                ek.append(np.sqrt(np.abs(self._reproject_slice(
                    var, cube.wcs, wcs_out, (ny, nx), flux_per_pixel))))
                wk.append(wl)
            order = np.argsort(wk)
            band_cubes.append({"name": name, "waves": np.asarray(wk)[order],
                               "flux": np.asarray(fk)[order],
                               "err": np.asarray(ek)[order]})
        band_cubes.sort(key=lambda b: float(b["waves"][0]))  # fix #2

        # 2. Flux-match bands onto a common scale (fix #1, prototype #3), by the
        #    chained walk or the joint rank-1 solve.
        if not self.flux_match:
            self.stitch_log = []
        elif self.flux_method == "joint":
            self.stitch_log = self._apply_joint_rank1(band_cubes)
        elif self.flux_method == "hybrid":
            self.stitch_log = self._apply_hybrid(band_cubes)
        else:
            self.stitch_log = self._apply_flux_matching(band_cubes)

        # 3. Concatenate, joining adjacent bands at their overlap midpoint so
        #    the noisy band roll-off edges are dropped (not kept-first).
        m = self._band_keep_masks(band_cubes)
        flux = np.concatenate([b["flux"][mk] for b, mk in zip(band_cubes, m)], axis=0)
        err = np.concatenate([b["err"][mk] for b, mk in zip(band_cubes, m)], axis=0)
        waves = np.concatenate([b["waves"][mk] for b, mk in zip(band_cubes, m)])
        bands = np.concatenate([np.full(int(mk.sum()), b["name"])
                                for b, mk in zip(band_cubes, m)])
        dq = (~np.isfinite(flux)).astype("float32")  # 1 where no coverage
        return StitchedCube(flux=flux, error=err, dq=dq, waves=waves,
                            wcs=wcs_out, bandnames=bands,
                            stitch_log=self.stitch_log,
                            err_rescale_log=self.err_rescale_log)

    @staticmethod
    def _bandname(cube):
        """A short label for a cube's sub-band, from its header if available."""
        hdr = getattr(cube, "general_header", None) or getattr(cube, "science_header", {})
        chan = hdr.get("CHANNEL") if hdr else None
        band = hdr.get("BAND") if hdr else None
        if chan is not None and band is not None:
            return f"ch{chan}_{band}"
        return getattr(cube, "filepath", "band").split("/")[-1].split(".")[0]
