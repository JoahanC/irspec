import matplotlib.pyplot as plt
import numpy as np 
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
from matplotlib.colors import Normalize, LogNorm
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from astropy.io import fits
from astropy.visualization.wcsaxes import add_beam, add_scalebar

from astropy.visualization import ZScaleInterval
from cubespec import CubeSpec
plt.rcParams["font.family"] = "sans-serif"

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_0_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB0", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name1 = "IR 23128-N 0"
spec_dict = {name1: {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}}
coords = {name1: (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])}

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_1_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB1", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name2 = "IR 23128-N 1"
spec_dict[name2] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name2] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_2_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB2", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name3 = "IR 23128-N 2"
spec_dict[name3] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name3] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])

scaling = [1, 50, 2000, 10000, 2e5, 5e6, 1e8]
plt.style.use('dark_background')
colors = ["orangered", "cyan", "magenta", "lime", "pink", "brown", "gold"]


hst_file1 = "./../archival_data/hubble_images/IR23128/j9cv82010_drz.fits"
hst_file2 = "./../archival_data/hubble_images/IR23128/j9cv82020_drz.fits"
miri_file1 = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"

from matplotlib.gridspec import GridSpec
plt.rcParams["font.family"] = "Verdana"

south_nuc_ra = "23h15m46.750s"
south_nuc_dec = "-59d03m15.60s"
sourthern_skycoord = SkyCoord(south_nuc_ra, south_nuc_dec, frame="icrs")
sourthern_nuc_rad = 1.5

fig = plt.figure()
fig.set_size_inches(12, 6)
gs1 = GridSpec(2,2)
gs1.update(left=-0.25, right=0.50, bottom=0, top=1, wspace=0.00,hspace=0.02)
names = ["IR23128a", "IR23128b", "IR23128c", "IR23128r0.8", "IR23128r1.1", "IR23128r1.5", "IR23128r1.8"]
linestyles = ["solid", "dashed", "dotted", "dashdot", "solid", "dashed", "dotted", "dashdot"]
linestyles = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid"]
files = [hst_file1]
wfc_pixscale = 0.05
for filename in files:
    hdu = fits.open(filename, ext=1)[1]
    wcs = WCS(hdu.header)
    ax1 = plt.subplot(gs1[0, :], projection=wcs)
    z = ZScaleInterval()
    z1, z2 = z.get_limits(hdu.data)
    #ax = plt.subplot()
    
    ax1.set_xlim(2000, 2300)
    ax1.set_ylim(2300, 2600)
    ax1.imshow(hdu.data, norm=LogNorm(0.01, 5), cmap="bone_r", origin="lower")
    
    for idx, name in enumerate(coords):
        coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
        x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
        rad = coords[name][2] / wfc_pixscale
        aperture = plt.Circle((x_pix, y_pix), rad, color=colors[idx], label=names[idx], fill=False)
        ax1.add_patch(aperture)
    
    
    #rad = sourthern_nuc_rad / wfc_pixscale
    #x_pix, y_pix = skycoord_to_pixel(sourthern_skycoord, wcs)
    #aperture = plt.Circle((x_pix, y_pix), rad, color="yellow", label="IR23128S", fill=False)
    #ax1.add_patch(aperture)
    
    #ax1.set_title()
    x_text_pos = 0.05 * (2300 - 2000) + 2000
    y_text_pos = 0.95 * (2600 - 2300) + 2300
    ax1.tick_params(axis='x', labelsize=0)
    ax1.tick_params(axis='y', labelsize=0)
    ax1.text(2005, 2570, "HST F435W", fontsize=20, color="black")
    ax1.set_xlabel("RA (J2000)", fontsize=0)
    ax1.set_ylabel("DEC (J2000)", fontsize=0)
    #plt.grid(True)


    #ax.set_xlim()

miri_names = ["a", "b", "c", "r0.8", "r1.1", "r1.5", "r1.8", "S"]
miri_xoffsets = [4, 4, 0, 5, 7, 9, 11, 3]
miri_yoffsets = [-1, 0, 4, 5, 7, 9, 11, 3]

miri_pixscale = 0.11 #asec per pix#
files = [miri_file1]
for filename in files:
    hdu = fits.open(filename, ext=1)[1]
    wcs = WCS(hdu.header)
    ax2 = plt.subplot(gs1[1, :], projection=wcs)
    z = ZScaleInterval()
    z1, z2 = z.get_limits(hdu.data)
    #ax = plt.subplot()
    #ax.set_xlim(1500, 2800)
    #ax.set_ylim(2000, 3000)
    z1 = 2
    z2 = 5000
    ax2.imshow(hdu.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    ax2.tick_params(axis='x', labelsize=0)
    ax2.tick_params(axis='y', labelsize=0)
    ax2.set_xlim(225, 325)
    ax2.set_ylim(225, 325)
    x_text_pos = 0.03 * (325 - 225) + 225
    y_text_pos = 0.90 * (325 - 225) + 225
    ax2.text(x_text_pos, y_text_pos, "JWST MIRI", fontsize=20, color="black")
    for idx, name in enumerate(coords):
        coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
        x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
        rad = coords[name][2] / miri_pixscale
        aperture = plt.Circle((x_pix, y_pix), rad, ls=linestyles[idx], color=colors[idx], label=names[idx], fill=False)
        ax2.add_patch(aperture)
        #ax2.text(x=x_pix + miri_xoffsets[idx],y=y_pix+miri_yoffsets[idx],s=miri_names[idx], color=colors[idx])
    
    #rad = sourthern_nuc_rad / miri_pixscale
    #x_pix, y_pix = skycoord_to_pixel(sourthern_skycoord, wcs)
    #aperture = plt.Circle((x_pix, y_pix), rad, color="yellow", label="IR23128S", fill=False)
    #ax2.add_patch(aperture)
    #ax2.text(x_pix + 10, y_pix + 10, s=miri_names[3], color="yellow")
    
    ax2.set_xlabel("RA (J2000)", fontsize=0)
    ax2.set_ylabel("DEC (J2000)", fontsize=0)
    #plt.grid(True)


#ax3 = plt.subplot(gs1[2, :],sharex=ax1)


gs2 = GridSpec(1, 1)
gs2.update(left=0.30, right=0.90, top=0.95, hspace=0.00)
ax4 = plt.subplot(gs2[0,0])
ax4.yaxis.set_label_position("right")
ax4.yaxis.tick_right()
for idx, key in enumerate(spec_dict):
    ax4.plot(spec_dict[key]["wave"], spec_dict[key]["flux"] * scaling[idx], linestyle = linestyles[idx], color=colors[idx], label=names[idx])
ax4.tick_params(direction='in', which='both', length=6, width=1, top=True)
ax4.tick_params(axis='x', labelsize=18)
ax4.tick_params(axis='y', labelsize=18)

ax4.set_xlabel(r"$\lambda_{rest}$ (micron)", fontsize=20)
ax4.set_ylabel(r'Scaled Flux Density', fontsize=18, rotation=270, labelpad=20)
#ax1.set_xscale('log')
ax4.set_yscale('log')
ax4.xaxis.set_minor_locator(MultipleLocator(1))
ax4.xaxis.grid(True, which='major', alpha=0.5)
ax4.xaxis.grid(True, which='minor', alpha=0.25)
ax4.legend(loc='best')
#ax4.set_aspect('equal', adjustable='box') 
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.savefig("ir23128n_spectra2.pdf", dpi=1000)




spec_obj = CubeSpec("./../", "param_files", "IR23128-N_0_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB0", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name1 = "IR 23128-N 0"
spec_dict = {name1: {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}}
coords = {name1: (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])}

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_1_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB1", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name2 = "IR 23128-N 1"
spec_dict[name2] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name2] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_2_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB2", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name3 = "IR 23128-N 2"
spec_dict[name3] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name3] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])

spec_obj = CubeSpec("./../", "param_files", "IR23128-N_12_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB12", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
name4 = "IR 23128-N 3"
spec_dict[name4] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
coords[name4] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])

scaling = [1, 50, 2000, 10000, 2e5, 5e6, 1e8]
plt.style.use('dark_background')
colors = ["orangered", "cyan", "magenta", "lime", "pink", "brown", "gold"]


hst_file1 = "./../archival_data/hubble_images/IR23128/j9cv82010_drz.fits"
hst_file2 = "./../archival_data/hubble_images/IR23128/j9cv82020_drz.fits"
miri_file1 = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"
nirspec_file1 = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"


from matplotlib.gridspec import GridSpec
plt.rcParams["font.family"] = "Verdana"

south_nuc_ra = "23h15m46.750s"
south_nuc_dec = "-59d03m15.60s"
sourthern_skycoord = SkyCoord(south_nuc_ra, south_nuc_dec, frame="icrs")
sourthern_nuc_rad = 1.5

fig = plt.figure()
fig.set_size_inches(12, 6)
gs1 = GridSpec(2,2)
gs1.update(left=-0.25, right=0.50, bottom=0, top=1, wspace=0.00,hspace=0.02)
names = ["IR23128a", "IR23128b", "IR23128c", "IR23128r1.8", "IR23128r1.1", "IR23128r1.5", "IR23128r1.8"]
linestyles = ["solid", "dashed", "dotted", "dashdot", "solid", "dashed", "dotted", "dashdot"]
linestyles = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid"]
files = [hst_file1]
wfc_pixscale = 0.05
for filename in files:
    hdu = fits.open(filename, ext=1)[1]
    wcs = WCS(hdu.header)
    ax1 = plt.subplot(gs1[0, :], projection=wcs)
    z = ZScaleInterval()
    z1, z2 = z.get_limits(hdu.data)
    #ax = plt.subplot()
    
    ax1.set_xlim(2000, 2300)
    ax1.set_ylim(2300, 2600)
    ax1.imshow(hdu.data, norm=LogNorm(0.01, 5), cmap="bone_r", origin="lower")
    
    for idx, name in enumerate(coords):
        coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
        x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
        rad = coords[name][2] / wfc_pixscale
        aperture = plt.Circle((x_pix, y_pix), rad, color=colors[idx], label=names[idx], fill=False)
        ax1.add_patch(aperture)
    
    
    #rad = sourthern_nuc_rad / wfc_pixscale
    #x_pix, y_pix = skycoord_to_pixel(sourthern_skycoord, wcs)
    #aperture = plt.Circle((x_pix, y_pix), rad, color="yellow", label="IR23128S", fill=False)
    #ax1.add_patch(aperture)
    
    #ax1.set_title()
    x_text_pos = 0.05 * (2300 - 2000) + 2000
    y_text_pos = 0.95 * (2600 - 2300) + 2300
    ax1.tick_params(axis='x', labelsize=0)
    ax1.tick_params(axis='y', labelsize=0)
    ax1.text(2005, 2570, "HST F435W", fontsize=20, color="black")
    ax1.set_xlabel("RA (J2000)", fontsize=0)
    ax1.set_ylabel("DEC (J2000)", fontsize=0)
    #plt.grid(True)


    #ax.set_xlim()

miri_names = ["a", "b", "c", "r0.8", "r1.1", "r1.5", "r1.8", "S"]
miri_xoffsets = [4, 4, 0, 5, 7, 9, 11, 3]
miri_yoffsets = [-1, 0, 4, 5, 7, 9, 11, 3]

miri_pixscale = 0.11 #asec per pix#
files = [miri_file1]
for filename in files:
    hdu = fits.open(filename, ext=1)[1]
    wcs = WCS(hdu.header)
    ax2 = plt.subplot(gs1[1, :], projection=wcs)
    z = ZScaleInterval()
    z1, z2 = z.get_limits(hdu.data)
    #ax = plt.subplot()
    #ax.set_xlim(1500, 2800)
    #ax.set_ylim(2000, 3000)
    z1 = 2
    z2 = 5000
    ax2.imshow(hdu.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    ax2.tick_params(axis='x', labelsize=0)
    ax2.tick_params(axis='y', labelsize=0)
    ax2.set_xlim(225, 325)
    ax2.set_ylim(225, 325)
    x_text_pos = 0.03 * (325 - 225) + 225
    y_text_pos = 0.90 * (325 - 225) + 225
    ax2.text(x_text_pos, y_text_pos, "JWST MIRI", fontsize=20, color="black")
    for idx, name in enumerate(coords):
        coord_obj = SkyCoord(coords[name][0], coords[name][1], unit="deg")
        x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
        rad = coords[name][2] / miri_pixscale
        aperture = plt.Circle((x_pix, y_pix), rad, ls=linestyles[idx], color=colors[idx], label=names[idx], fill=False)
        ax2.add_patch(aperture)
        #ax2.text(x=x_pix + miri_xoffsets[idx],y=y_pix+miri_yoffsets[idx],s=miri_names[idx], color=colors[idx])
    
    #rad = sourthern_nuc_rad / miri_pixscale
    #x_pix, y_pix = skycoord_to_pixel(sourthern_skycoord, wcs)
    #aperture = plt.Circle((x_pix, y_pix), rad, color="yellow", label="IR23128S", fill=False)
    #ax2.add_patch(aperture)
    #ax2.text(x_pix + 10, y_pix + 10, s=miri_names[3], color="yellow")
    
    ax2.set_xlabel("RA (J2000)", fontsize=0)
    ax2.set_ylabel("DEC (J2000)", fontsize=0)
    #plt.grid(True)


#ax3 = plt.subplot(gs1[2, :],sharex=ax1)


gs2 = GridSpec(1, 1)
gs2.update(left=0.30, right=0.90, top=0.95, hspace=0.00)
ax4 = plt.subplot(gs2[0,0])
ax4.yaxis.set_label_position("right")
ax4.yaxis.tick_right()
for idx, key in enumerate(spec_dict):
    ax4.plot(spec_dict[key]["wave"], spec_dict[key]["flux"] * scaling[idx], linestyle = linestyles[idx], color=colors[idx], label=names[idx])
ax4.tick_params(direction='in', which='both', length=6, width=1, top=True)
ax4.tick_params(axis='x', labelsize=18)
ax4.tick_params(axis='y', labelsize=18)

ax4.set_xlabel(r"$\lambda_{rest}$ (micron)", fontsize=20)
ax4.set_ylabel(r'Scaled Flux Density', fontsize=18, rotation=270, labelpad=20)
#ax1.set_xscale('log')
ax4.set_yscale('log')
ax4.xaxis.set_minor_locator(MultipleLocator(1))
ax4.xaxis.grid(True, which='major', alpha=0.5)
ax4.xaxis.grid(True, which='minor', alpha=0.25)
ax4.legend(loc='best')
#ax4.set_aspect('equal', adjustable='box') 
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.savefig("ir23128n_spectra3.pdf", dpi=1000)