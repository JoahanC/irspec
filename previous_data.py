import numpy as np 
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import fits
from matplotlib.colors import Normalize, LogNorm
from cubespec import CubeSpec
from astropy.visualization import ZScaleInterval

plt.style.use('dark_background')
spec_obj = CubeSpec("./../", "param_files", "IR23128-N_0_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB0", mode="SB")
#north_cubes = spec_obj.perform_single_extraction_custom()
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name1 = "IR 23128-N 0"
spec_dict = {name1: {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}}
print(spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"])
coords = {name1: ("23h15m46.689s", "-59d03m10.80s", spec_obj.param_dict["user_r_ap"])}

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_1_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB1", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name2 = "IR 23128-N 1"
spec_dict[name2] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
print(spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"])

coords[name2] = ("23h15m46.725s", "-59d03m10.20s", spec_obj.param_dict["user_r_ap"])

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_2_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB2", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name3 = "IR 23128-N 2"
spec_dict[name3] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
print(spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"])

coords[name3] = ("23h15m46.785s", "-59d03m09.87s", spec_obj.param_dict["user_r_ap"])

"""spec_obj = CubeSpec("./../", "param_files", "IR23128-N_12_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB12", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name4 = "IR 23128-N 3"
spec_dict[name4] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name4] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])

spec_obj = CubeSpec("./../", "param_files", "IR23128-S_6_single_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN6", mode="AGN")
#south_cubes = spec_obj.perform_single_extraction_custom()
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name5 = "IR 23128-S 4"
spec_dict[name5] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name5] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])"""

hst_file1 = "./../archival_data/hubble_images/IR23128/j9cv82010_drz.fits"
hst_file2 = "./../archival_data/hubble_images/IR23128/j9cv82020_drz.fits"
miri_file1 = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"

names = ["IR23128a", "IR23128b", "IR23128c", "IR23128r0.8", "IR23128r1.1", "IR23128r1.5", "IR23128r1.8"]
linestyles = ["solid", "dashed", "dotted", "dashdot", "solid", "dashed", "dotted", "dashdot"]
linestyles = ["solid", "solid", "solid", "dashdot", "solid", "dashed", "dotted", "dashdot"]
colors = ["orangered", "cyan", "magenta", "lime", "pink", "brown", "gold"]

miri_pixscale = 0.11
hdu = fits.open(miri_file1, ext=1)[1]
wcs = WCS(hdu.header)
fig = plt.figure()
fig.set_size_inches(6, 6)
ax = plt.subplot(projection=wcs)
z1 = 5
z2 = 3000
z = ZScaleInterval()
ax.imshow(hdu.data, norm=LogNorm(z1, z2), origin='lower', cmap="plasma")
for idx, name in enumerate(coords):
    coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
    x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
    rad = coords[name][2] / miri_pixscale
    aperture = plt.Circle((x_pix, y_pix), rad, lw=2, ls=linestyles[idx], color=colors[idx], label=names[idx], fill=False)
    ax.add_patch(aperture)
ax.set_xlim(215, 335)
ax.set_ylim(215, 335)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.set_xlabel("RA (J2000)", fontsize=20)
ax.set_ylabel("Dec (J2000)", fontsize=20)
ax.set_title("MIRI F770W", loc="right", fontsize=22)
plt.grid(alpha=0.25)
plt.savefig("miri_444_spots.pdf", dpi=1000, bbox_inches='tight')
plt.close()

spitzer_3_file = "./../archival_data/spitzer_images/SPITZER_I1_12313344_0000_7_E8700443_maic.fits"
spitzer_8_file = "./../archival_data/spitzer_images/SPITZER_I4_12313344_0000_7_E8700562_maic.fits"

hdu = fits.open(spitzer_3_file, ext=1)
header = fits.getheader(spitzer_3_file)
wcs = WCS(header)
fig = plt.figure()
fig.set_size_inches(6, 6)
ax = plt.subplot(projection=wcs)
z1 = 100
z2 = 3000
z = ZScaleInterval()

data = fits.getdata(spitzer_3_file)[420:520,1100:1200]
#data = fits.getdata(spitzer_3_file)[520:420,1200:1100]

z1, z2 = z.get_limits(data)
print(z1, z2)
ax.imshow(data, norm=LogNorm(0.2, 30), cmap="plasma")
#for idx, name in enumerate(coords):
#    coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
#    x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
#    rad = coords[name][2] / miri_pixscale
#    aperture = plt.Circle((x_pix, y_pix), rad, lw=2, ls=linestyles[idx], color=colors[idx], label=names[idx], fill=False)
#    ax.add_patch(aperture)
print(np.shape(data))
#ax.set_xlim(1100, 1200)
#ax.set_ylim(420, 520)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.set_xlabel("RA (J2000)", fontsize=20)
ax.set_ylabel("Dec (J2000)", fontsize=20)
ax.set_title(r"IRAC 3.6$\mu m$", loc="right", fontsize=22)
plt.grid(alpha=0.25)
plt.savefig("spitzer_36.pdf", dpi=1000, bbox_inches="tight")
plt.close()
#plt.savefig("miri_444_spots.pdf", dpi=1000)

"""hdu = fits.open(spitzer_8_file, ext=1)
header = fits.getheader(spitzer_8_file)
wcs = WCS(header)
fig = plt.figure()
fig.set_size_inches(6, 6)
ax = plt.subplot(projection=wcs)
z1 = 100
z2 = 3000
z = ZScaleInterval()

data = fits.getdata(spitzer_8_file)#[420:520,1100:1200]
z1, z2 = z.get_limits(data)
print(z1, z2)
ax.imshow(data, norm=LogNorm(0.2, 200), origin='lower', cmap="bone")
#for idx, name in enumerate(coords):
#    coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
#    x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
#    rad = coords[name][2] / miri_pixscale
#    aperture = plt.Circle((x_pix, y_pix), rad, lw=2, ls=linestyles[idx], color=colors[idx], label=names[idx], fill=False)
#    ax.add_patch(aperture)
print(np.shape(data))
#ax.set_xlim(1100, 1200)
#ax.set_ylim(420, 520)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.set_xlabel("RA (J2000)", fontsize=20)
ax.set_ylabel("Dec (J2000)", fontsize=20)
ax.set_title(r"IRAC 8.0$\mu m$", loc="right", fontsize=22)
plt.grid(alpha=0.25)
plt.show()
#plt.savefit("spitzer_8-.pdf", dpi=1000, bbox_inches="tight")
#plt.savefig("miri_444_spots.pdf", dpi=1000)

wfc_pixescale=0.05"""