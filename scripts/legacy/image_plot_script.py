import glob
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import Normalize, LogNorm
from matplotlib.gridspec import GridSpec
from irspec.datacube import Datacube
from irspec.cubespec import CubeSpec
from regions import PolygonSkyRegion, PixCoord, CircleSkyRegion
from astropy.visualization import ZScaleInterval
from astropy.visualization.wcsaxes import add_beam, add_scalebar
from astropy import units as u
import matplotlib.font_manager as fm
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel

from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp


image_plots = False 
extraction_plots = True

test_path = glob.glob(f"./irspec/data/ch4-long_s3d.fits")[0]
testcube = Datacube(test_path, redshift=0.044601, verbose=True)

g395h_path = glob.glob(f"./irspec/data/g395h-f290lp_s3d.fits")[0]
g395hcube = Datacube(g395h_path, redshift=0.044601, verbose=True)
g395hregion_data = g395hcube.science_header["S_REGION"].split()
g395h_pixel_x = [6, -3, 82, 92]
g395h_pixel_y = [-3, 71, 83, 10]
g395h_pixcoord = PixCoord(g395h_pixel_x, g395h_pixel_y)
g395h_skycoord = g395h_pixcoord.to_sky(g395hcube.wcs)

ch1s_path = glob.glob(f"./irspec/data/ch1-short_s3d.fits")[0]
ch1scube = Datacube(ch1s_path, redshift=0.044601, verbose=True)
ch1sregion_data = ch1scube.science_header["S_REGION"].split()
ch1s_pixel_x = [26, 0, 22, 48]
ch1s_pixel_y = [0, 28, 48, 21]
ch1s_pixcoord = PixCoord(ch1s_pixel_x, ch1s_pixel_y)
ch1s_skycoord = ch1s_pixcoord.to_sky(ch1scube.wcs)

ch4l_path = glob.glob(f"./irspec/data/ch4-long_s3d.fits")[0]
ch4lcube = Datacube(ch4l_path, redshift=0.044601, verbose=True)
ch4region_data = ch4lcube.science_header["S_REGION"].split()
ch4l_pixel_x = [15, 35, 19, -1]
ch4l_pixel_y = [35, 15, -1, 19]
ch4l_pixcoord = PixCoord(ch4l_pixel_x, ch4l_pixel_y)
ch4l_skycoord = ch4l_pixcoord.to_sky(ch4lcube.wcs)

nch1s_path = glob.glob(f"./irspec/north_data/ch1-short_s3d.fits")[0]
nch1scube = Datacube(nch1s_path, redshift=0.044601, verbose=True)
nch1sregion_data = nch1scube.science_header["S_REGION"].split()
nch1s_pixel_x = [26, 0, 22, 48]
nch1s_pixel_y = [0, 28, 48, 21]
nch1s_pixcoord = PixCoord(nch1s_pixel_x, nch1s_pixel_y)
nch1s_skycoord = nch1s_pixcoord.to_sky(nch1scube.wcs)

nch4l_path = glob.glob(f"./irspec/north_data/ch4-long_s3d.fits")[0]
nch4lcube = Datacube(nch4l_path, redshift=0.044601, verbose=True)
nch4region_data = nch4lcube.science_header["S_REGION"].split()
nch4l_pixel_x = [15, 35, 19, -1]
nch4l_pixel_y = [35, 15, -1, 19]
nch4l_pixcoord = PixCoord(nch4l_pixel_x, nch4l_pixel_y)
nch4l_skycoord = nch4l_pixcoord.to_sky(nch4lcube.wcs)

hdu1 = fits.open(get_pkg_data_filename('galactic_center/gc_2mass_k.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('galactic_center/gc_msx_e.fits'))[0]



nircam_image_path = "./irspec/image_data/nircam/IR23128_f200w.fits"
miri_image_path = "./irspec/image_data/miri/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"
hst_image_path = "./irspec/image_data/hubble/j9cv82020_drz.fits"

south_nuc_ra = "23h15m46.750s"
south_nuc_dec = "-59d03m15.60s"
sourthern_skycoord = SkyCoord(south_nuc_ra, south_nuc_dec, frame="icrs")
sourthern_nuc_rad = 1.5

colors_n = ["blue", "magenta", "lime", "pink", "brown", "gold"]
colors_s = ["green", "orange", "red"]
scaling_n = [0.1, 10, 100, 1000, 10000, 1e5]
scaling_s = [30, 800, 100, 1e5]
# NIRSPEC
#scaling_s = [10, 30, 100, 1e5]
# PAH
#scaling_n = [0.1, 2, 100, 1000, 10000, 1e5]
#scaling_s = [1, 4, 10, 1e5]

if image_plots:

    # Plots individual image of NIRCAM + SOUTH IFS
    fig = plt.figure()
    hdu = fits.open(nircam_image_path, ext=1)[1]
    wcs = WCS(hdu.header)
    #ch1s_skycoord_plot = ch1s_skycoord.to_pixel(wcs)
    g395h_sky_region = PolygonSkyRegion(g395h_skycoord)
    g395h_pixel_region = g395h_sky_region.to_pixel(wcs)
    ch1s_sky_region = PolygonSkyRegion(ch1s_skycoord)
    ch1s_pixel_region = ch1s_sky_region.to_pixel(wcs)
    ch4l_sky_region = PolygonSkyRegion(ch4l_skycoord)
    ch4l_pixel_region = ch4l_sky_region.to_pixel(wcs)
    #ax2 = plt.subplot(projection=wcs)
    ax = plt.subplot()
    z = ZScaleInterval()
    z1, z2 = z.get_limits(hdu.data)
    z1 = 2
    z2 = 50
    ax.imshow(hdu.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    ch1s_pixel_region.plot(ax=ax, color="green", lw=2.0, label="CH1-short")
    ch4l_pixel_region.plot(ax=ax, color="red", lw=2.0, label="CH4-long")
    g395h_pixel_region.plot(ax=ax, color="blue", lw=2.0, label="g395h-f290")
    ax.tick_params(axis='x', labelsize=0)
    ax.tick_params(axis='y', labelsize=0)
    ax.set_xlim(3400, 3800)
    ax.set_ylim(3000, 3400)
    x_text_pos = 0.03 * (350 - 150) + 150
    y_text_pos = 0.90 * (400 - 200) + 200
    ax.text(x_text_pos, y_text_pos, "NIRCam F200W", fontsize=20, color="black")
    ax.set_xlabel("RA (J2000)", fontsize=0)
    ax.set_ylabel("DEC (J2000)", fontsize=0)
    ax.legend()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.savefig("./../../images/nircam200_ifs_s.png", dpi=300)
    plt.close()

    # Plots individual image of NIRCAM + NORTH IFS
    fig = plt.figure()
    hdu = fits.open(nircam_image_path, ext=1)[1]
    wcs = WCS(hdu.header)
    #ch1s_skycoord_plot = ch1s_skycoord.to_pixel(wcs)
    nch1s_sky_region = PolygonSkyRegion(nch1s_skycoord)
    nch1s_pixel_region = nch1s_sky_region.to_pixel(wcs)
    nch4l_sky_region = PolygonSkyRegion(nch4l_skycoord)
    nch4l_pixel_region = nch4l_sky_region.to_pixel(wcs)
    #ax2 = plt.subplot(projection=wcs)
    ax = plt.subplot()
    z = ZScaleInterval()
    z1, z2 = z.get_limits(hdu.data)
    z1 = 2
    z2 = 50
    ax.imshow(hdu.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    nch1s_pixel_region.plot(ax=ax, color="green", lw=2.0, label="CH1-short")
    nch4l_pixel_region.plot(ax=ax, color="red", lw=2.0, label="CH4-long")
    ax.tick_params(axis='x', labelsize=0)
    ax.tick_params(axis='y', labelsize=0)
    ax.set_xlim(3400, 3800)
    ax.set_ylim(3000, 3400)
    x_text_pos = 0.03 * (350 - 150) + 150
    y_text_pos = 0.90 * (400 - 200) + 200
    ax.text(x_text_pos, y_text_pos, "NIRCam F200W", fontsize=20, color="black")
    ax.set_xlabel("RA (J2000)", fontsize=0)
    ax.set_ylabel("DEC (J2000)", fontsize=0)
    ax.legend()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.savefig("./../../images/nircam200_ifs_n.png", dpi=300)
    plt.close()

    hdu_miri = fits.open(miri_image_path, ext=1)[1]
    hdu_nircam = fits.open(nircam_image_path, ext=1)[1]
    hdu_hst = fits.open(hst_image_path, ext=1)[1]
    wcs_miri = WCS(hdu_miri.header)
    wcs_nircam = WCS(hdu_nircam.header)
    wcs_hst = WCS(hdu_hst.header)

    # Plots side by side of MIRI and NIRCAM
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    gs1 = GridSpec(ncols=2,nrows=1)
    ax1 = fig.add_subplot(gs1[0], projection=wcs_miri)
    ax2 = fig.add_subplot(gs1[1], projection=wcs_nircam)
    z1 = 2
    z2 = 5000
    ax1.imshow(hdu_miri.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    z1 = 2
    z2 = 50
    ax2.imshow(hdu_nircam.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    #print(np.shape(hdu.data))
    ax1.set_xlabel("Right Ascension (J2000)")
    ax2.set_xlabel("Right Ascension (J2000)")
    ax1.set_ylabel("Declination (J2000)")
    ax2.yaxis.set_visible(False)
    ax1.set_xlim(50, 510)
    ax1.set_ylim(50, 510)
    ax2.set_xlim(3000, 4200)
    ax2.set_ylim(2600, 3800)
    ax1.text(310, 310, "N", fontsize=20)
    ax1.text(230, 230, "S", fontsize=20)
    ax2.text(3400, 3345, "N", fontsize=20)
    ax2.text(3686, 3100, "S", fontsize=20)
    ax1.set_title("JWST MIRI F770W")
    ax2.set_title("JWST NIRCam F200W")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 10 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=20, family='Helvetica')
    add_scalebar(ax1, scalebar_angle, label="10 kpc", color="black", fontproperties=fontprops)
    add_scalebar(ax2, scalebar_angle, label="10 kpc", color="black", fontproperties=fontprops)
    ax1.grid()
    ax2.grid()
    plt.savefig("./../../images/miri_ifs.png", dpi=300)
    plt.close()

    # Plots side by side of HST and MIRI
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    gs1 = GridSpec(ncols=2,nrows=1)
    ax1 = fig.add_subplot(gs1[0])
    ax2 = fig.add_subplot(gs1[1])

    z1 = 0.01
    z2 = 2000
    # [y:y][x:x]
    # [1250:2750] [200:1500,200:1500]
    ax1.imshow(np.flipud(np.fliplr(hdu_hst.data[1500:2800,1380:2680][200:1500,200:1500])), norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    z1 = 1
    z2 = 3000
    ax2.imshow(hdu_miri.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    #print(np.shape(hdu.data))
    #ax1.set_xlabel("Right Ascension (J2000)")
    #ax1.set_ylabel("Declination (J2000)")
    ax2.set_xlim(50, 510)
    ax2.set_ylim(50, 510)
    #ax2.set_xlim(3000, 4200)
    #ax2.set_ylim(2600, 3800)
    ax1.text(600, 500, "N", fontsize=20)
    ax1.text(400, 250, "S", fontsize=20)
    ax2.text(310, 310, "N", fontsize=20)
    ax2.text(230, 230, "S", fontsize=20)
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 10 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=20, family='Helvetica')
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    #add_scalebar(ax1, scalebar_angle, label="10 kpc", color="black", fontproperties=fontprops)
    #add_scalebar(ax2, scalebar_angle, label="10 kpc", color="black", fontproperties=fontprops)
    #ax1.grid()
    #ax2.grid()
    ax1.set_title("HST WFC F770W")
    ax2.set_title("JWST MIRI F770W")
    plt.savefig("./../../images/hst_miri.png", dpi=300)
    plt.close()

    # Plots cutout of NIRCAM nuclei
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    ax = fig.add_subplot(projection=wcs_nircam)
    z1 = 2
    z2 = 50
    ax.imshow(hdu_nircam.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    #print(np.shape(hdu.data))
    ax.set_xlabel("Right Ascension (J2000)")
    ax.set_ylabel("Declination (J2000)")
    ax.set_xlim(3400, 3800)
    ax.set_ylim(3000, 3400)
    ax.text(3500, 3345, "N", fontsize=28)
    ax.text(3686, 3100, "S", fontsize=28)
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 1 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=20, family='Helvetica')
    add_scalebar(ax, scalebar_angle, label="1 kpc", color="black", fontproperties=fontprops)
    ax.grid()
    plt.savefig("./../../images/nircam_zoom.png", dpi=300)
    plt.close()

    # NIRCAM + IFS EXTRACTIONS SIDE BY SIDE
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    gs1 = GridSpec(ncols=2,nrows=1)
    ax1 = fig.add_subplot(gs1[0], projection=wcs_nircam)
    ax2 = fig.add_subplot(gs1[1], projection=wcs_nircam)
    nch1s_sky_region = PolygonSkyRegion(nch1s_skycoord)
    nch1s_pixel_region = nch1s_sky_region.to_pixel(wcs)
    nch4l_sky_region = PolygonSkyRegion(nch4l_skycoord)
    nch4l_pixel_region = nch4l_sky_region.to_pixel(wcs)
    ch1s_sky_region = PolygonSkyRegion(ch1s_skycoord)
    ch1s_pixel_region = ch1s_sky_region.to_pixel(wcs)
    ch4l_sky_region = PolygonSkyRegion(ch4l_skycoord)
    ch4l_pixel_region = ch4l_sky_region.to_pixel(wcs)

    z1 = 2
    z2 = 50
    ax1.imshow(hdu_nircam.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    z1 = 2
    z2 = 50
    ax2.imshow(hdu_nircam.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
    nch1s_pixel_region.plot(ax=ax1, color="green", lw=2.0, label="CH1-short")
    nch4l_pixel_region.plot(ax=ax1, color="red", lw=2.0, label="CH4-long")
    ch1s_pixel_region.plot(ax=ax2, color="green", lw=2.0, label="CH1-short")
    ch4l_pixel_region.plot(ax=ax2, color="red", lw=2.0, label="CH4-long")

    ax1.set_xlabel("Right Ascension (J2000)")
    ax2.set_xlabel("Right Ascension (J2000)")
    ax1.set_ylabel("Declination (J2000)")
    ax2.yaxis.set_visible(False)
    ax1.set_xlim(3400, 3800)
    ax1.set_ylim(3000, 3400)
    ax2.set_xlim(3400, 3800)
    ax2.set_ylim(3000, 3400)
    ax1.set_title("Northern IFS")
    ax2.set_title("Southern IFS")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 1 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=20, family='Helvetica')
    add_scalebar(ax1, scalebar_angle, label="1 kpc", color="black", fontproperties=fontprops)
    add_scalebar(ax2, scalebar_angle, label="1 kpc", color="black", fontproperties=fontprops)
    ax1.grid()
    ax2.grid()
    plt.savefig("./../../images/sidebyside_ifs.png", dpi=300)
    plt.close()

if extraction_plots:
    # Load all extractions
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                    "ir23128s_single1S_param.txt", "./../../creta_extractions/single1S",
                    "./../../cafe_output/single1S_AGN", 0.044601, testcube, mode="AGN")
    af = spec_obj.open_asdf()
    wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
    flux = np.asarray(af['cafefit']['obsspec']['flux'])
    flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
    name1 = "Sa"
    spec_dict_s = {name1: {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}}
    coords_s = {name1: (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])}
    
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files", 
                    "ir23128s_single3S_param.txt", "./../../creta_extractions/single3S",
                    "./../../cafe_output/single3S_AGN", 0.044601, testcube, mode="AGN")
    af = spec_obj.open_asdf()
    wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
    flux = np.asarray(af['cafefit']['obsspec']['flux'])
    flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
    name2 = "Sb"
    spec_dict_s[name2] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
    coords_s[name2] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])
    
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files", 
                    "ir23128s_single2S_param.txt", "./../../creta_extractions/single2S",
                    "./../../cafe_output/single2S_SB", 0.044601, testcube, mode="AGN")
    af = spec_obj.open_asdf()
    wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
    flux = np.asarray(af['cafefit']['obsspec']['flux'])
    flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
    name3 = "Sc"
    spec_dict_s[name3] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
    coords_s[name3] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])
    
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/north_data/", "./irspec/param_files/", 
                    "ir23128s_single1N_param.txt", "./../../creta_extractions/single1N",
                    "./../../cafe_output/single1N_SB", 0.044601, testcube, mode="SB")
    af = spec_obj.open_asdf()
    wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
    flux = np.asarray(af['cafefit']['obsspec']['flux'])
    flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
    name1 = "Na"
    spec_dict_n = {name1: {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}}
    coords_n = {name1: (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])}
    
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/north_data/", "./irspec/param_files/", 
                    "ir23128s_single2N_param.txt", "./../../creta_extractions/single2N",
                    "./../../cafe_output/single2N_SB", 0.044601, testcube, mode="SB")
    af = spec_obj.open_asdf()
    wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
    flux = np.asarray(af['cafefit']['obsspec']['flux'])
    flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
    name2 = "Nb"
    spec_dict_n[name2] = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
    coords_n[name2] = (spec_obj.param_dict["user_ra"], spec_obj.param_dict["user_dec"], spec_obj.param_dict["user_r_ap"])
    
    
    
    
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    gs1 = GridSpec(2,2)
    gs1.update(left=-0.25, right=0.50, bottom=0, top=1, wspace=0.00,hspace=0.02)
    names_n = ["Na", "Nb", "Nc", "IR23128Sr1.25", "IR23128Sr1.5", "IR23128Sr1.7"]
    names_s = ["Sa", "Sb", "Sc", "IR23128Sr1.25", "IR23128Sr1.5", "IR23128Sr1.7"]
    #colors_s = []
    #colors_n = []
    linestyles = ["solid", "dashed", "dotted", "dashdot", "solid", "dashed", "dotted", "dashdot"]
    linestyles = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid"]
    radii_n = [1.0, 0.5, 0.5]
    radii_s = [1.2, 0.5, 0.5]
    
    
    
    
    
    files = [nircam_image_path]
    wfc_pixscale = 0.05
    for filename in files:
        hdu = fits.open(filename, ext=1)[1]
        wcs = WCS(hdu.header)
        
        nch1s_sky_region = PolygonSkyRegion(nch1s_skycoord)
        nch1s_pixel_region = nch1s_sky_region.to_pixel(wcs)
        nch4l_sky_region = PolygonSkyRegion(nch4l_skycoord)
        nch4l_pixel_region = nch4l_sky_region.to_pixel(wcs)
        
        
        ax1 = plt.subplot(gs1[0, :], projection=wcs)
        z = ZScaleInterval()
        #z1, z2 = z.get_limits(hdu.data)
        z1 = 2
        z2 = 50
        #ax = plt.subplot()
        
        #ax1.set_xlim(2000, 2300)
        #ax1.set_ylim(2300, 2600)
        ax1.imshow(hdu.data, norm=LogNorm(z1, z2), cmap="bone_r", origin="lower")
        
        for idx, name in enumerate(coords_n):
            coord_obj = SkyCoord(coords_n[name][0], coords_n[name][1], unit="deg")
            aper_region = CircleSkyRegion(coord_obj, radii_n[idx]*u.arcsecond)
            aper_pixel_region = aper_region.to_pixel(wcs)
            aper_pixel_region.plot(ax=ax1, color=colors_n[idx], lw=2.0, label="CH1-short")
            
        #for idx, name in enumerate(coords_n):
        #    coord_obj = SkyCoord(coords_n[name][0], coords_n[name][1], unit="deg")
        #    x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
        #    rad = coords_n[name][2] / wfc_pixscale
        #    aperture = plt.Circle((x_pix, y_pix), rad, color=colors_n[idx], label=names_n[idx], fill=False)
        #    ax1.add_patch(aperture)
        
        
        #rad = sourthern_nuc_rad / wfc_pixscale
        #x_pix, y_pix = skycoord_to_pixel(sourthern_skycoord, wcs)
        #aperture = plt.Circle((x_pix, y_pix), rad, color="yellow", label="IR23128S", fill=False)
        #ax1.add_patch(aperture)
        nch1s_pixel_region.plot(ax=ax1, color="green", lw=2.0, label="CH1-short")
        #nch4l_pixel_region.plot(ax=ax1, color="red", lw=2.0, label="CH4-long")
        #ax1.set_title()
        #x_text_pos = 0.05 * (2300 - 2000) + 2000
        #y_text_pos = 0.95 * (2600 - 2300) + 2300
        #ax1.tick_params(axis='x', labelsize=0)
        #ax1.tick_params(axis='y', labelsize=0)
        ax1.set_xlim(3440, 3640)
        ax1.set_ylim(3180, 3380)
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        #ax1.text(2005, 2570, "HST F435W", fontsize=20, color="black")
        ax1.set_xlabel("RA (J2000)", fontsize=0)
        ax1.set_ylabel("DEC (J2000)", fontsize=0)
    
    miri_names = ["r0.5", "r0.75", "r1.0", "r1.25", "r1.5", "r1.7"]
    miri_xoffsets = [4, 4, 0, 5, 7, 9, 11, 3]
    miri_yoffsets = [-1, 0, 4, 5, 7, 9, 11, 3]

    miri_pixscale = 0.11 #asec per pix#
    files = [nircam_image_path]
    for filename in files:
        hdu = fits.open(filename, ext=1)[1]
        wcs = WCS(hdu.header)
        
        g395h_sky_region = PolygonSkyRegion(g395h_skycoord)
        g395h_pixel_region = g395h_sky_region.to_pixel(wcs)
        ch1s_sky_region = PolygonSkyRegion(ch1s_skycoord)
        ch1s_pixel_region = ch1s_sky_region.to_pixel(wcs)
        ch4l_sky_region = PolygonSkyRegion(ch4l_skycoord)
        ch4l_pixel_region = ch4l_sky_region.to_pixel(wcs)
        
        ax2 = plt.subplot(gs1[1, :], projection=wcs)
        z = ZScaleInterval()
        z1, z2 = z.get_limits(hdu.data)
        #ax = plt.subplot()
        #ax.set_xlim(1500, 2800)
        #ax.set_ylim(2000, 3000)
        z1 = 2
        z2 = 50
        ax2.imshow(hdu.data, norm=LogNorm(z1, z2), origin='lower', cmap="bone_r")
        #ax2.tick_params(axis='x', labelsize=0)
        #ax2.tick_params(axis='y', labelsize=0)
        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        #ax2.set_xlim(225, 325)
        #ax2.set_ylim(225, 325)
        ax2.set_xlim(3500, 3700)
        ax2.set_ylim(3050, 3250)
        ch1s_pixel_region.plot(ax=ax2, color="green", lw=2.0, label="CH1-short")
        #ch4l_pixel_region.plot(ax=ax2, color="red", lw=2.0, label="CH4-long")
        g395h_pixel_region.plot(ax=ax2, color="orange", lw=2.0, label="g395h-f290")
        #x_text_pos = 0.03 * (325 - 225) + 225
        #y_text_pos = 0.90 * (325 - 225) + 225
        #ax2.text(x_text_pos, y_text_pos, "JWST MIRI", fontsize=20, color="black")
        for idx, name in enumerate(coords_s):
            coord_obj = SkyCoord(coords_s[name][0], coords_s[name][1], unit="deg")
            aper_region = CircleSkyRegion(coord_obj, radii_s[idx]*u.arcsecond)
            aper_pixel_region = aper_region.to_pixel(wcs)
            aper_pixel_region.plot(ax=ax2, color=colors_s[idx], lw=2.0, label="CH1-short")
            #x_pix, y_pix = skycoord_to_pixel(coord_obj, wcs)
            #rad = coords[name][2] / miri_pixscale
            #aperture = plt.Circle((x_pix, y_pix), rad, ls=linestyles[idx], color=colors[idx], label=names[idx], fill=False)
            #ax2.add_patch(aperture)
            #ax2.text(x=x_pix + miri_xoffsets[idx],y=y_pix+miri_yoffsets[idx],s=miri_names[idx], color=colors[idx])
        

        
        ax2.set_xlabel("RA (J2000)", fontsize=0)
        ax2.set_ylabel("DEC (J2000)", fontsize=0)
    
    
    gs2 = GridSpec(1, 1)
    gs2.update(left=0.30, right=0.90, top=0.95, hspace=0.00)
    ax4 = plt.subplot(gs2[0,0])
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    for idx, key in enumerate(spec_dict_n):
        ax4.plot(spec_dict_n[key]["wave"], spec_dict_n[key]["flux"] * scaling_n[idx], linestyle = linestyles[idx], color=colors_n[idx], label=names_n[idx])
    for idx, key in enumerate(spec_dict_s):
        ax4.plot(spec_dict_s[key]["wave"], spec_dict_s[key]["flux"] * scaling_s[idx], linestyle = linestyles[idx], color=colors_s[idx], label=names_s[idx])
    ax4.tick_params(direction='in', which='both', length=6, width=1, top=True)
    ax4.tick_params(axis='x', labelsize=18)
    ax4.tick_params(axis='y', labelsize=18)

    ax4.set_ylim(0.001, 800)
    # NIRSPEC FULL
    #ax4.set_xlim(2.85, 5)
    #ax4.set_ylim(0.05, 40)
    # NIRSPEC
    #ax4.set_xlim(4.5, 5)
    #ax4.set_ylim(0.05, 0.8)
    # PAHS
    #ax4.set_xlim(6, 13)
    #ax4.set_ylim(0.001, 1)
    # LATE STAGE
    ax4.set_xlim(13, 27)
    ax4.set_ylim(0.001, 1000)
    #ax4.axvspan(6.1, 6.5, color="gray", alpha=0.3, label="PAH")
    #ax4.axvspan(7.2, 8.1, color="gray", alpha=0.3)
    #ax4.axvspan(8.4, 8.9, color="gray", alpha=0.3)
    #ax4.axvspan(11.0, 11.6, color="gray", alpha=0.3)
    #ax4.axvspan(12.5, 12.9, color="gray", alpha=0.3)
    ax4.legend(loc='upper center', prop={'size': 12})
    
    ax4.set_xlabel(r"$\lambda_{rest}$ (micron)", fontsize=20)
    ax4.set_ylabel(r'Scaled Flux Density', fontsize=18, rotation=270, labelpad=20)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    #ax4.xaxis.set_minor_locator(MultipleLocator(1))
    ax4.xaxis.grid(True, which='major', alpha=0.5)
    ax4.xaxis.grid(True, which='minor', alpha=0.25)
    ax4.legend(loc='best', prop={'size': 12})
    #ax4.set_aspect('equal', adjustable='box') 
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.savefig("ir23128s_spectra_latestage.png", dpi=1000)