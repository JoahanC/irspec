"""
This file defines a datacube object which stores any metadata and real
data associated with a JWST stage 3 datacube product.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.constants as const


class Datacube:
    """An object which stores metadata of a Level3 processed JWST datacube.

    Attributes:
        filename (str): The location of the JWST datacube.
        general_header (Header): The non-science high level header of the
            datacube.
        science_header (Header): The header containing the observation
            parameters.
        science_data (ndarray): The flux array for IFS measurements.
        error_data (ndarray): The flux error array for IFS measurements.
        dq_data (ndarray): The data quality flag array for IFS measurements.
    """


    def __init__(self, filepath, redshift=0, verbose=False):
        """Inits Datacube.

        Args:
            filepath (str): The path to the datacube.
            redshift (float, optional): The redshift of the observations
                contained in the datacube.
            verbose (bool, optional): Whether to display logs.
        """
        self.filepath = filepath
        self.redshift = redshift
        self._verbose = verbose
        self._read_fits()


    ### Housekeeping Methods


    @staticmethod
    def _as_quantity(value, unit):
        """Coerces a plain number to a Quantity in `unit`; passes Quantities through."""
        if isinstance(value, u.Quantity):
            return value
        return value * u.Unit(unit)


    def _read_fits(self):
        """Reads and stores header and data information from the FITS file."""
        with fits.open(self.filepath) as hdul:
            self.general_header = hdul[0].header.copy()
            self.science_header = hdul[1].header.copy()
            self.science_data = hdul[1].data.copy()
            self.error_data = hdul[2].data.copy()
            self.dq_data = hdul[3].data.copy()

        self.im_shape = np.shape(self.science_data[0])
        self.ax3_len = len(self.science_data)

        self.wcs = WCS(self.science_header).celestial
        self.plot_header = self.wcs.to_header()

        # Housekeeping variables for computations
        self._wv_unit = u.Unit(self.science_header["CUNIT3"])
        self._sb_unit = u.Unit(self.science_header["BUNIT"])
        self._ref_wv = self.science_header["CRVAL3"] * self._wv_unit
        self._delta_wv = self.science_header["CDELT3"] * self._wv_unit
        self.wvs = (self._ref_wv + self._delta_wv * np.arange(self.ax3_len)) / (1 + self.redshift)
        self.area_sr = self.science_header["PIXAR_SR"] * u.steradian

    # Datacube inspection methods

    def spaxel_values(self, x_pix, y_pix):
        """Returns the spaxel surface brightness spectrum, spectrum errors,
        and spectrum data quality flag at a given (x, y) coordinate.

        Note: if BUNIT is "Jy/pix" rather than a surface-brightness unit
        (e.g. "MJy/sr"), the returned quantities are converted to an
        integrated flux (Jy) instead, since "per pixel" already implies
        a fixed aperture equal to one spaxel.
        """
        sb = self.science_data[:, x_pix, y_pix] * self._sb_unit
        sb_err = self.error_data[:, x_pix, y_pix] * self._sb_unit
        if self._sb_unit == u.Unit("Jy/pix"):
            sb = (sb * u.pixel).to(u.Jy)
            sb_err = (sb_err * u.pixel).to(u.Jy)
        return sb, sb_err, self.dq_data[:, x_pix, y_pix]

    def mrs_resolving_power(self, wavelength):
        """Calculate the MRS resolving power at a given wavelength"""
        if isinstance(wavelength, u.Quantity):
            return 4603 - 128 * wavelength.value
        return 4603 - 128 * wavelength

    # Datacube computation methods

    def wv_to_idx(self, wv, unit="um"):
        """Returns the closest index corresponding to a given wavelength value.

        Args:
            wv (float or Quantity): The wavelength to locate.
            unit (Unit): The unit of ``wv`` if it is not already a Quantity.

        Raises:
            ValueError: If the wavelength is out of bounds.
        """
        wv = self._as_quantity(wv, unit).to(self.wvs[0].unit)
        if wv >= self.wvs[0] and wv < self.wvs[-1]:
            return np.argmin(abs(self.wvs - wv))
        else:
            raise ValueError("Wavelength is out of bounds!")

    def wv_to_vel(self, wv, ref_wv, unit="um"):
        """Returns the velocity shift corresponding to a wavelength shift.

        Args:
            wv (float or Quantity): The wavelength at which to evaluate the
                velocity shift.
            ref_wv (float or Quantity): The reference wavelength to evaluate
                the wavelength shift.
            unit (Unit): The units of the provided wavelengths.

        Returns:
            Quantity: The velocity shift in kilometers/second.
        """
        wv = self._as_quantity(wv, unit)
        ref_wv = self._as_quantity(ref_wv, unit)
        return (const.c * (wv - ref_wv) / ref_wv).to(u.kilometer / u.second)

    def vel_to_wv(self, vel, ref_wv, vel_unit="km/s", wv_unit="um"):
        """Returns the wavelength corresponding to a velocity shift and a
        reference wavelength.

        Args:
            vel (float or Quantity): The velocity shift in kilometers/second.
            ref_wv (float or Quantity): The reference wavelength to evaluate
                the velocity shift.
            vel_unit (Unit): The units of the velocity shift.
            wv_unit (Unit): The units of the reference wavelength and of the
                returned wavelength.

        Returns:
            Quantity: The wavelength shift in the specified ``wv_unit``.
        """
        ref_wv = self._as_quantity(ref_wv, wv_unit)
        vel = self._as_quantity(vel, vel_unit)
        return (ref_wv * vel / const.c).to(u.micron)

    def vel_to_sigma(self, vel_disp, ref_wavelength, vel_unit="km/s", wv_unit="um"):
        """Converts a velocity-space dispersion (standard deviation) to a
        wavelength-space dispersion via the Doppler relation.

        Args:
            vel_disp (float or Quantity): The velocity-space dispersion
                (standard deviation).
            ref_wavelength (float or Quantity): The reference wavelength.
            vel_unit (Unit): The units of ``vel_disp``, if not already a
                Quantity.
            wv_unit (Unit): The units of ``ref_wavelength``, if not already a
                Quantity.

        Returns:
            Quantity: The wavelength-space dispersion, in microns.
        """
        vel_disp = self._as_quantity(vel_disp, vel_unit)
        ref_wavelength = self._as_quantity(ref_wavelength, wv_unit)
        return (ref_wavelength * vel_disp / const.c).to(u.micron)

    def sigma_to_disp(self, wv_sigma, line_center, unit="um"):
        """Converts a wavelength-space dispersion (standard deviation) to a
        velocity-space dispersion via the Doppler relation.

        Note:
            This converts between standard deviations (sigma), not FWHMs. To
            start from a FWHM, first divide by 2*sqrt(2*ln(2)) (~2.3548) to
            get sigma.

        Args:
            wv_sigma (float or Quantity): The wavelength-space dispersion
                (standard deviation).
            line_center (float or Quantity): The reference wavelength.
            unit (Unit): The units of ``wv_sigma`` and ``line_center``, if not
                already Quantities.

        Returns:
            Quantity: The velocity-space dispersion, in kilometers/second.
        """
        wv_sigma = self._as_quantity(wv_sigma, unit)
        line_center = self._as_quantity(line_center, unit)
        return (const.c * (wv_sigma / line_center)).to(u.kilometer / u.second)

    ### Visualization Methods

    def display_spaxel_spectra(self, x_pix, y_pix):
        """Displays the de-redshifted spectrum of a spaxel."""
        sb, sb_err, dq = self.spaxel_values(x_pix, y_pix)
        fig, ax = plt.subplots()
        ax.errorbar(self.wvs, sb, sb_err, capsize=2, color="black", ecolor="black")
        ax.set_xlabel("Wavelength [μm]")
        ax.set_ylabel("Surface Brightness [MJy/sr]")
        ax.set_title(f"({x_pix}, {y_pix})", loc="right")
        plt.show()

    def _plot_image(self, data, cbar_label, cmap="plasma", scale=None):
        """Plots a 2D image array with a labeled colorbar. Returns fig, ax."""
        fig, ax = plt.subplots()
        image = ax.imshow(data, norm=scale, origin="lower", cmap=cmap)
        cax = fig.colorbar(image)
        cax.set_label(cbar_label, fontsize=14, rotation=270, labelpad=20)
        ax.set_xlabel("XPIX")
        ax.set_ylabel("YPIX")
        return fig, ax

    def display_flux_slice(self, idx, scale="log"):
        fig, ax = self._plot_image(self.science_data[idx], "Surface Brightness [MJy/sr]", scale=scale)
        ax.set_title(f"Index: {idx}", loc="left")
        ax.set_title(f"Wavelength: {round(self.wvs[idx].value, 3)} μm", loc="right")
        plt.show()

    def display_error_slice(self, idx, scale="log"):
        self._plot_image(self.error_data[idx], "Surface Brightness Error [MJy/sr]", scale=scale)
        plt.show()

    def display_dq_slice(self, idx):
        """Displays the data quality array for a given array slice."""
        self._plot_image(self.dq_data[idx], "Data Quality Flags", cmap="Set1")
        plt.show()

    def display_dq_complete(self):
        """Displays each spaxel's percentage of good data quality."""
        qualities = np.zeros(np.shape(self.dq_data[0]))
        for slice in self.dq_data:
            qualities += slice == 0
        qualities /= len(self.dq_data) / 100

        self._plot_image(qualities, "Percentage of non-null data", cmap="binary")
        plt.show()
