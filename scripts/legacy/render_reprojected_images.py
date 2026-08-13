import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle
import os as _os, sys as _sys
# plotparams lives in the irspec package under irspec/src; this script runs
# from irspec/, so put the package root on the path before importing it.
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "src"))
from irspec.plotparams import PlotParams

import astropy.units as u
import numpy as np
from astropy.wcs import *
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.visualization import ZScaleInterval
from astropy.visualization.wcsaxes import add_scalebar
from astropy.stats import sigma_clipped_stats

plt.rcParams["font.family"] = "Helvetica"
z = ZScaleInterval()


hdu1 = fits.open(get_pkg_data_filename('galactic_center/gc_2mass_k.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('galactic_center/gc_msx_e.fits'))[0]

#hst_file1 = "./../archival_data/hubble_images/IR23128/j9cv82010_drz.fits"
hst_file1 = "./src/irspec/image_data/hubble/j9cv82010_drz.fits"

hst_file2 = "./../archival_data/hubble_images/IR23128/j9cv82020_drz.fits"
#miri_file1 = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"
miri_file1 = "./src/irspec/image_data/miri/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"

#nircam_file1 = "./../archival_data/nircam_images/IR23128/IR23128_f200w.fits"
nircam_file1 = "./src/irspec/image_data/nircam/IR23128_f200w.fits"

nircam_file2 = "./../archival_data/nircam_images/IR23128/IR23128_f356w.fits"
#spitzer_3_file = "./../archival_data/spitzer_images/SPITZER_I1_12313344_0000_7_E8700443_maic.fits"

# ALMA band 6 aggregate continuum, 247.6 GHz (1.21mm), project 2017.1.00759.S.
# 0.298" x 0.239" beam; the 40" field covers both the north and south pointings.
alma_cont_file = "./src/irspec/image_data/alma/member.uid___A001_X1289_X55c.ESO148-IG002_sci.spw25_27_29_31.cont.I.pbcor.fits"
alma_pb_file = "./src/irspec/image_data/alma/member.uid___A001_X1289_X55c.ESO148-IG002_sci.spw25_27_29_31.cont.I.pb.fits.gz"

# Distances to IR 23128. Every scalebar below wants the angular-diameter
# distance, not the luminosity distance: D_A = D_L/(1+z)^2 is the Etherington
# reciprocity relation, so deriving one from the other keeps these bars tied
# to whatever cosmology produced the 194.99 Mpc rather than assuming one.
#
# D_A is the smaller of the two (178.7 Mpc here), so a given physical length
# subtends a *larger* angle than the luminosity distance implies. Drawing
# these bars from D_L made them about 9% too short for their labels, and any
# size read off one came out overestimated by that much.
REDSHIFT = 0.044601
LUMINOSITY_DISTANCE = 194.99 * u.Mpc
ANGULAR_DIAMETER_DISTANCE = LUMINOSITY_DISTANCE / (1 + REDSHIFT) ** 2

datacube_north_1 = f"./src/irspec/north_data/ch1-short_s3d.fits"
datacube_north_4 = f"./src/irspec/north_data/ch4-long_s3d.fits"
datacube_south_1 = f"./src/irspec/data/ch1-short_s3d.fits"
datacube_south_4 = f"./src/irspec/data/ch4-long_s3d.fits"

"""
hdu = fits.open(spitzer_3_file, ext=1)
data = fits.getdata(spitzer_3_file)

#array, footprint = reproject_interp(hdu2, hdu1.header)
wcs_out, shape_out = find_optimal_celestial_wcs(hdu)

plt.subplot(projection=wcs_out)
plt.imshow(data)#, origin='lower')
plt.grid(color='white', ls='solid')
plt.xlabel('ra')
plt.ylabel('dec')
plt.show()"""
plt.style.use("default")


def crop_flip(array, x_low, x_high, y_low, y_high):
    """Crop and flip an array reprojected onto hdu_ref, the way every panel does."""
    return np.fliplr(np.flipud(array[y_low:y_high, x_low:x_high]))


def load_alma_continuum(cont_file, pb_file):
    """Load the ALMA continuum map as a 2D HDU and measure its noise level.

    The pipeline products carry four axes (RA, Dec, FREQ, STOKES), which
    reproject cannot take, so the degenerate axes are dropped and the WCS
    reduced to its celestial part. This delivery ships only the primary-beam
    corrected image, whose noise rises toward the field edge, so the flat-noise
    map used for the RMS has to be rebuilt as pbcor * pb.
    """
    hdu = fits.open(cont_file)[0]
    image = np.squeeze(hdu.data)
    pb = np.squeeze(fits.open(pb_file)[0].data)

    flat_noise = image * pb
    _, _, rms = sigma_clipped_stats(flat_noise[np.isfinite(flat_noise)],
                                    sigma=3.0, maxiters=5)

    # Past the half-power point the primary-beam correction amplifies noise
    # faster than signal, so those pixels only contribute spurious contours.
    image = np.where(pb > 0.5, image, np.nan)

    header = WCS(hdu.header).celestial.to_header()
    for key in ("BMAJ", "BMIN", "BPA", "BUNIT"):
        header[key] = hdu.header[key]
    return fits.PrimaryHDU(image, header), rms


def add_alma_contours(ax, array, rms, x_low, x_high, y_low, y_high,
                      sigmas=(5, 10, 20, 40, 80), color="cyan"):
    """Overlay ALMA continuum contours on a panel already drawn with imshow.

    The contours go through the same crop and flip as the background image, so
    they register with what is displayed regardless of the panel's WCS.

    The lowest level is 5 sigma rather than the usual 3: the map holds some
    2e4 independent beams, so a 3 sigma floor litters the field with a few
    dozen noise islands.
    """
    ax.contour(crop_flip(array, x_low, x_high, y_low, y_high),
               levels=np.array(sigmas) * rms, colors=color,
               linewidths=0.8, alpha=0.9)




hdu_ref = fits.open(hst_file1)[1]
hdu_miri = fits.open(miri_file1)[1]
hdu_nircam = fits.open(nircam_file1, ext=1)[1]
#hdu_spitzer = fits.open(spitzer_3_file)[0]

hdu_n_ch1 = fits.open(datacube_north_1)[1]
hdu_s_ch1 = fits.open(datacube_south_1)[1]

hdu_n_ch1_new_header = hdu_n_ch1.header
hdu_s_ch1_new_header = hdu_s_ch1.header

hdu_n_ch1_new_header["NAXIS"] = 2
hdu_n_ch1_new_header['WCSAXES'] = 2
del hdu_n_ch1_new_header['NAXIS3']
del hdu_n_ch1_new_header['CRPIX3']
del hdu_n_ch1_new_header['CRVAL3']
del hdu_n_ch1_new_header['CTYPE3']
del hdu_n_ch1_new_header['CUNIT3']
del hdu_n_ch1_new_header['CDELT3']
del hdu_n_ch1_new_header['PC1_3']
del hdu_n_ch1_new_header['PC2_3']
del hdu_n_ch1_new_header['PC3_1']
del hdu_n_ch1_new_header['PC3_2']
del hdu_n_ch1_new_header['PC3_3']

hdu_s_ch1_new_header["NAXIS"] = 2
hdu_s_ch1_new_header['WCSAXES'] = 2
del hdu_s_ch1_new_header['NAXIS3']
del hdu_s_ch1_new_header['CRPIX3']
del hdu_s_ch1_new_header['CRVAL3']
del hdu_s_ch1_new_header['CTYPE3']
del hdu_s_ch1_new_header['CUNIT3']
del hdu_s_ch1_new_header['CDELT3']
del hdu_s_ch1_new_header['PC1_3']
del hdu_s_ch1_new_header['PC2_3']
del hdu_s_ch1_new_header['PC3_1']
del hdu_s_ch1_new_header['PC3_2']
del hdu_s_ch1_new_header['PC3_3']

data_arr = hdu_n_ch1.data[100]
data_arr2 = hdu_s_ch1.data[100]

fits.writeto('output_file.fits', data_arr, hdu_n_ch1_new_header, overwrite=True)
fits.writeto('output_file2.fits', data_arr2, hdu_s_ch1_new_header, overwrite=True)


hdu_test = fits.open("output_file.fits")
hdu_test2 = fits.open("output_file2.fits")
#print(repr(hdu_n_ch1_new_header))

mrs_n_ch1_array, mrs_n_ch1_footprint = reproject_interp(hdu_test, hdu_ref.header)
mrs_s_ch1_array, mrs_s_ch1_footprint = reproject_interp(hdu_test2, hdu_ref.header)


gc_distance = ANGULAR_DIAMETER_DISTANCE
scalebar_length = 10 * u.kpc
scalebar_angle = (scalebar_length / gc_distance).to(
    u.deg, equivalencies=u.dimensionless_angles()
)
fontprops = fm.FontProperties(size=18, family='Helvetica')

x_low_spitzer = 1550
x_high_spitzer = 2650
y_low_spitzer = 1800
y_high_spitzer = 2900

x_low_miri = 1950
x_high_miri = 2350
y_low_miri = 2200
y_high_miri = 2600

x_low_mrs_n = 2075
x_high_mrs_n = 2225
y_low_mrs_n = 2315
y_high_mrs_n = 2465

x_low_mrs_s = 2100
x_high_mrs_s = 2250
y_low_mrs_s = 2400
y_high_mrs_s = 2550

miri_len = x_high_miri - x_low_miri
mrs_len = 150




miri_array, miri_footprint = reproject_interp(hdu_miri, hdu_ref.header)
nircam_array, nircam_footprint = reproject_interp(hdu_nircam, hdu_ref.header)

hdu_alma, alma_rms = load_alma_continuum(alma_cont_file, alma_pb_file)
alma_array, alma_footprint = reproject_interp(hdu_alma, hdu_ref.header)
print(f"ALMA continuum RMS: {alma_rms * 1e6:.1f} uJy/beam")

"""
spitzer_array, spitzer_footprint = reproject_interp(hdu_spitzer, hdu_ref.header)



"""

"""rec = Rectangle((y_high_spitzer - y_high_miri, x_high_spitzer - x_high_miri), miri_len, miri_len, fill=False, color="gold", lw=2, label="MIRI")


fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.subplot(projection=WCS(hdu_ref.header))
ax1.imshow(np.fliplr(np.flipud(spitzer_array[y_low_spitzer:y_high_spitzer, x_low_spitzer:x_high_spitzer])), origin='lower', norm=LogNorm(0.1, 100), cmap="plasma")
ax1.coords['ra'].set_axislabel('RA (J2000)', fontsize=22)
ax1.coords['dec'].set_axislabel('Dec (J2000)', fontsize=22)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title(r'Spitzer IRAC 3.6$\mu$m', loc="right", fontsize=28)
#ax1.add_patch(rec)
#ax1.legend(fontsize=20)
add_scalebar(ax1, scalebar_angle, label="10 kpc", color="white", fontproperties=fontprops)
plt.savefig("./irac_3_final_miri_less.pdf", dpi=1200, bbox_inches="tight")
plt.close()

rec = Rectangle((y_high_spitzer - y_high_miri, x_high_spitzer - x_high_miri), miri_len, miri_len, fill=False, color="gold", lw=2, label="MIRI")

fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.subplot(projection=WCS(hdu_ref.header))
ax1.imshow(np.fliplr(np.flipud(hdu_ref.data[y_low_spitzer:y_high_spitzer, x_low_spitzer:x_high_spitzer])), origin='lower', norm=LogNorm(0.01, 20), cmap="plasma")
ax1.coords['ra'].set_axislabel('RA (J2000)', fontsize=22)
ax1.coords['dec'].set_axislabel('Dec (J2000)', fontsize=22)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title(r'Hubble WFC F814W', loc="right", fontsize=28)
add_scalebar(ax1, scalebar_angle, label="10 kpc", color="white", fontproperties=fontprops)
#ax1.add_patch(rec)
#ax1.legend(fontsize=20)
plt.savefig("./hst_final_miri_less.pdf", dpi=1200, bbox_inches="tight")
plt.close()"""
ppram = PlotParams(palatte="light", scaling="presentation", grid=True)
plt.rcParams["font.family"] = "Helvetica"

recn = Rectangle((y_high_miri - y_high_mrs_n + 60, x_high_miri - x_high_mrs_n + 20), mrs_len*0.5, mrs_len*0.8, fill=False, angle=30, color="gold", lw=2, label="MRS North")
recs = Rectangle((y_high_miri - y_high_mrs_s + 120, x_high_miri - x_high_mrs_s - 40), mrs_len*0.5, mrs_len*0.8, fill=False, angle=30, color="lime", lw=2, label="MRS South")


fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.subplot(projection=WCS(hdu_ref.header))
ax1.imshow(np.fliplr(np.flipud(miri_array[y_low_miri:y_high_miri, x_low_miri:x_high_miri])), origin='lower', norm=LogNorm(6, 3000), cmap="plasma")
ax1.coords['ra'].set_axislabel('RA (J2000)', fontsize=22)
ax1.coords['dec'].set_axislabel('Dec (J2000)', fontsize=22)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title('MIRI F770W', loc="right", fontsize=28)
ax1.add_patch(recn)
ax1.add_patch(recs)
add_alma_contours(ax1, alma_array, alma_rms,
                  x_low_miri, x_high_miri, y_low_miri, y_high_miri)
ax1.plot([], [], color="cyan", lw=1.5, label="ALMA 1.2mm")
ax1.legend(fontsize=20)
add_scalebar(ax1, scalebar_angle, label="10 kpc", color="white", fontproperties=fontprops)#
#plt.show()
plt.savefig("./miri_f700.png", dpi=1200, bbox_inches="tight")

gc_distance = ANGULAR_DIAMETER_DISTANCE
scalebar_length = 1 * u.kpc
scalebar_angle = (scalebar_length / gc_distance).to(
    u.deg, equivalencies=u.dimensionless_angles()
)

fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.subplot(projection=WCS(hdu_ref.header))
ax1.imshow(np.fliplr(np.flipud(mrs_n_ch1_array[y_low_mrs_n:y_high_mrs_n, x_low_mrs_n:x_high_mrs_n])), origin='lower', norm=LogNorm(20, 200), cmap="plasma")
ax1.coords['ra'].set_axislabel('RA (J2000)', fontsize=22)
ax1.coords['dec'].set_axislabel('Dec (J2000)', fontsize=22)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title('MRS CH1 North', loc="right", fontsize=28)
add_alma_contours(ax1, alma_array, alma_rms,
                  x_low_mrs_n, x_high_mrs_n, y_low_mrs_n, y_high_mrs_n)
#ax1.add_patch(rec)
#ax1.legend(fontsize=20)
# Foreground, not white: the MRS cubes cover only part of this panel and the
# rest is bare figure background, so a white bar drawn here sits on white.
add_scalebar(ax1, scalebar_angle, label="1 kpc", color=ppram.foreground, fontproperties=fontprops)
plt.savefig("./mrs_north.png", dpi=1200, bbox_inches="tight")
plt.close()

fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.subplot(projection=WCS(hdu_ref.header))
ax1.imshow(np.fliplr(np.flipud(mrs_s_ch1_array[y_low_mrs_s:y_high_mrs_s, x_low_mrs_s:x_high_mrs_s])), origin='lower', norm=LogNorm(20, 200), cmap="plasma")
ax1.coords['ra'].set_axislabel('RA (J2000)', fontsize=22)
ax1.coords['dec'].set_axislabel('Dec (J2000)', fontsize=22)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_title('MRS CH1 South', loc="right", fontsize=28)
add_alma_contours(ax1, alma_array, alma_rms,
                  x_low_mrs_s, x_high_mrs_s, y_low_mrs_s, y_high_mrs_s)
#ax1.add_patch(rec)
#ax1.legend(fontsize=20)
# Foreground, not white: the MRS cubes cover only part of this panel and the
# rest is bare figure background, so a white bar drawn here sits on white.
add_scalebar(ax1, scalebar_angle, label="1 kpc", color=ppram.foreground, fontproperties=fontprops)
#plt.show()
plt.savefig("./mrs_south.png", dpi=1200, bbox_inches="tight")
plt.close()

x_low_nircam = 2050
x_high_nircam = 2250
y_low_nircam = 2325
y_high_nircam = 2525

x_low_nircam = 2000
x_high_nircam = 2300
y_low_nircam = 2275
y_high_nircam = 2575

gc_distance = ANGULAR_DIAMETER_DISTANCE
scalebar_length = 10 * u.kpc
scalebar_angle = (scalebar_length / gc_distance).to(
    u.deg, equivalencies=u.dimensionless_angles()
)

fig = plt.figure()
fig.set_size_inches(10, 10)
ax = plt.subplot(projection=WCS(hdu_ref.header))
ax.imshow(np.fliplr(np.flipud(nircam_array[y_low_nircam:y_high_nircam, x_low_nircam:x_high_nircam])), origin='lower', norm=LogNorm(0.5, 1000), cmap="plasma")
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlabel("RA (J2000)", fontsize=22)
ax.set_ylabel("Dec (J2000)", fontsize=22)
ax.set_title("NIRCAM F200W", loc="right", fontsize=28)
add_alma_contours(ax, alma_array, alma_rms,
                  x_low_nircam, x_high_nircam, y_low_nircam, y_high_nircam)
ax.plot([], [], color="cyan", lw=1.5, label="ALMA 1.2mm")
ax.legend(fontsize=20)
add_scalebar(ax, scalebar_angle, label="10 kpc", color="white", fontproperties=fontprops)#

plt.savefig("nircam_f200w.png", dpi=1000, bbox_inches='tight')