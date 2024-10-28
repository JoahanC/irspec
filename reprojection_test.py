""" 
This file contains an example for reprojecting HST/Spitzer/JWST images.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle

import astropy.units as u
import numpy as np
from astropy.wcs import *
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.visualization import ZScaleInterval
from astropy.visualization.wcsaxes import add_scalebar

plt.rcParams["font.family"] = "Helvetica"
z = ZScaleInterval()


hdu1 = fits.open(get_pkg_data_filename('galactic_center/gc_2mass_k.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('galactic_center/gc_msx_e.fits'))[0]

hst_file1 = "./../archival_data/hubble_images/IR23128/j9cv82010_drz.fits"
hst_file2 = "./../archival_data/hubble_images/IR23128/j9cv82020_drz.fits"
miri_file1 = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"
nircam_file1 = "./../archival_data/nircam_images/IR23128/IR23128_f200w.fits"
nircam_file2 = "./../archival_data/nircam_images/IR23128/IR23128_f356w.fits"
spitzer_3_file = "./../archival_data/spitzer_images/SPITZER_I1_12313344_0000_7_E8700443_maic.fits"

datacube_north_1 = f"./../input_data/IR23128-N/iras23128-n_ch1-short_s3d.fits"
datacube_north_4 = f"./../input_data/IR23128-N/iras23128-n_ch4-long_s3d.fits"
datacube_south_1 = f"./../input_data/IR23128-S/Level3_ch1-short_s3d.fits"
datacube_south_4 = f"./../input_data/IR23128-S/Level3_ch4-long_s3d.fits"

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
plt.style.use("dark_background")

hdu_ref = fits.open(hst_file1)[1]
hdu_miri = fits.open(miri_file1)[1]
hdu_spitzer = fits.open(spitzer_3_file)[0]

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

"""mrs_n_ch1_array, mrs_n_ch1_footprint = reproject_interp(hdu_test, hdu_ref.header)
mrs_s_ch1_array, mrs_s_ch1_footprint = reproject_interp(hdu_test2, hdu_ref.header)"""


gc_distance = 194.99 * u.Mpc
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
ax1.set_title('JWST MIRI F770W', loc="right", fontsize=28)
ax1.add_patch(recn)
ax1.add_patch(recs)
ax1.legend(fontsize=20)
add_scalebar(ax1, scalebar_angle, label="10 kpc", color="white", fontproperties=fontprops)#
#plt.show()
plt.savefig("./miri_final_mrs_on.pdf", dpi=1200, bbox_inches="tight")

"""gc_distance = 194.99 * u.Mpc
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
ax1.set_title(r'JWST MIRI MRS North', loc="right", fontsize=28)
#ax1.add_patch(rec)
#ax1.legend(fontsize=20)
add_scalebar(ax1, scalebar_angle, label="1 kpc", color="white", fontproperties=fontprops)
plt.savefig("./mrs_north.pdf", dpi=1200, bbox_inches="tight")
plt.close()

fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.subplot(projection=WCS(hdu_ref.header))
ax1.imshow(np.fliplr(np.flipud(mrs_s_ch1_array[y_low_mrs_s:y_high_mrs_s, x_low_mrs_s:x_high_mrs_s])), origin='lower', norm=LogNorm(20, 200), cmap="plasma")
ax1.coords['ra'].set_axislabel('RA (J2000)', fontsize=22)
ax1.coords['dec'].set_axislabel('Dec (J2000)', fontsize=22)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title(r'JWST MIRI MRS South', loc="right", fontsize=28)
#ax1.add_patch(rec)
#ax1.legend(fontsize=20)
add_scalebar(ax1, scalebar_angle, label="1 kpc", color="white", fontproperties=fontprops)
#plt.show()
plt.savefig("./mrs_south.pdf", dpi=1200, bbox_inches="tight")
plt.close()"""