
import numpy as np 
import os
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
from astropy.io import ascii
from astropy.visualization import ZScaleInterval
from plotparams import PlotParams
from cubespec import CubeSpec
from astropy.visualization.wcsaxes import add_scalebar
import matplotlib.font_manager as fm
import scipy.special as sp

def gaussian_integral(amplitude, center, sigma, x_1, x_2):
    return amplitude * (sp.erf((center - x_1)) / (np.sqrt(2) * sigma) - sp.erf((center - x_2)) / (np.sqrt(2) * sigma)) / 2

def render_multicomponent_plot(data, wcs, savefig=False):
    """
    Renders a spaxel map illustrating the number of components fit to 
    each spaxel.
    """
    
    # Initialize base aarray
    base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))
    for idx, _ in enumerate(data["XPIX"]):
        if np.isnan(data["G3AMP"][idx]):
            if np.isnan(data["G2AMP"][idx]):
                if np.isnan(data["G1AMP"][idx]):
                    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 0
                else:
                    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 1
            else:
                base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 2
        else:
            base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 3
    
    fig = plt.figure()
    ax = plt.subplot(projection=wcs)
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('bone', np.max(base_array) - np.min(base_array) + 1)
    image = ax.imshow(base_array, cmap=cmap, vmin=np.min(base_array) - 0.5, vmax=np.max(base_array) + 0.5, origin="lower")
    ax.scatter(21, 24, c="white", edgecolors="black", marker="*", s=1000)
    plt.colorbar(image, ticks=np.arange(np.min(base_array), np.max(base_array) + 1))
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title("# of Components", loc="right")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 5 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
    if savefig:
        plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/compnum.png", dpi=600)
        plt.close()
    else:
        plt.show()


def render_totflux_plot(data, wcs, savefig=False):
    """
    Renders a spaxel map illustrating the number of components fit to 
    each spaxel
    """
    
    # Initialize base aarray
    base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))
    
    for idx, _ in enumerate(data["XPIX"]):
        flux_val = 0
        if not np.isnan(data["G1AMP"][idx]):
            
            g1_gamma = np.abs(data["G1SIGMA"][idx]) * 2.355 / data["G1CEN"][idx]
            g1_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (data["G1AMP"][idx] * u.Jy * g1_gamma / (data["G1CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
            flux_val += g1_flux.value
            
            if not np.isnan(data["G2AMP"][idx]):
                g2_gamma = np.abs(data["G2SIGMA"][idx]) * 2.355 / data["G2CEN"][idx]
                g2_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (data["G2AMP"][idx] * u.Jy * g2_gamma / (data["G2CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                flux_val += g2_flux.value
                if not np.isnan(data["G3AMP"][idx]):
                    g3_gamma = np.abs(data["G3SIGMA"][idx]) * 2.355 / data["G3CEN"][idx]
                    g3_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (data["G3AMP"][idx] * u.Jy * g3_gamma / (data["G3CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                    flux_val += g3_flux.value
                    if flux_val < 0:
                        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 0
                    else:
                        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = flux_val
                else:

                    if flux_val < 0:
                        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 0
                    else:
                        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = flux_val
            else:

                if flux_val < 0:
                    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = 0
                else:
                    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = flux_val
        else:
            base_array[data["XPIX"][idx]][data["YPIX"][idx]] = flux_val

    fig = plt.figure()
    ax = plt.subplot(projection=wcs)
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    image = ax.imshow(base_array, cmap=cmap, norm=LogNorm(), origin="lower")
    ax.scatter(21, 24, c="white", edgecolors="black", marker="*", s=1000)
    cax = plt.colorbar(image)
    cax.set_label(r"[W/m$^2$]", fontsize=24, rotation=270, labelpad=25)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title("Total Flux", loc="right")
    ax.set_title(line_name, loc="left")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 5 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
    if savefig:
        plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/flux.png", dpi=600)
        plt.close()
    else:
        plt.show()


def render_amplitude_plot(data, line_name, wcs, param="G1AMP", savefig=False):
    """
    Renders a spaxel map illustrating the number of components fit to 
    each spaxel
    """
    
    # Initialize base array
    base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))
    fluxes = []
    for idx, _ in enumerate(data["XPIX"]):
        pref_idx = param[1]
        gamma = np.abs(data["G"+pref_idx+"SIGMA"][idx]) * 2.355 / data["G"+pref_idx+"CEN"][idx]
        flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (data[param][idx] * u.Jy * gamma / (data["G"+pref_idx+"CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = flux.value
        fluxes.append(flux.value)
    
    max_flux = np.nanpercentile(fluxes, 99)
    
    if "1" in param:
        comp_name = r"$F_{1}$"
    if "2" in param:
        comp_name = r"$F_{2}$"
    if "3" in param:
        comp_name = r"$F_{3}$"
    
    fig = plt.figure()
    ax = plt.subplot(projection=wcs)
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    image = ax.imshow(base_array, cmap=cmap, norm=LogNorm(vmax=max_flux), origin="lower")
    ax.scatter(21, 24, c="white", edgecolors="black", marker="*", s=1000)
    cax = plt.colorbar(image)
    cax.set_label(r"[W/m$^2$]", fontsize=24, rotation=270, labelpad=25)

    #ax.set_facecolor("gray")
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title(comp_name, loc="right")
    ax.set_title(line_name, loc="left")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 5 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
    if savefig:
        plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/{param}.png", dpi=600)
        plt.close()
    else:
        plt.show()
        
def render_rel_vel_plot(data, line_center, line_name, wcs, param="G2CEN", savefig=False):
    """
    Renders a spaxel map illustrating the number of components fit to 
    each spaxel
    """
    
    # Initialize base aarray
    max_val = 0
    min_val = 0
    base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))
    all_relvels = []
    for idx, _ in enumerate(data["XPIX"]):
        #print(line_center, data[param][idx])
        
        rel_vel = ((const.c * (data[param][idx] - line_center)/line_center).to(u.kilometer / u.second)).value
        if rel_vel > max_val:
            max_val = rel_vel
        if rel_vel < min_val:
            min_val = rel_vel
        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = rel_vel
        all_relvels.append(rel_vel)

    
    high_percentile = np.nanpercentile(all_relvels, 97.5)
    
    if "1" in param:
        comp_name = r"$v_{1}$"
    if "2" in param:
        comp_name = r"$v_{2}$"
    if "3" in param:
        comp_name = r"$v_{3}$"
    
    max_bound = np.max([max_val, np.abs(min_val)])
    max_bound = high_percentile
    max_bound = 800
    fig = plt.figure()
    ax = plt.subplot(projection=wcs)
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('RdBu_r')
    image = ax.imshow(base_array, vmin=-max_bound, vmax=max_bound, cmap=cmap, origin="lower")
    ax.scatter(21, 24, c="white", edgecolors="black", marker="*", s=1000)
    cax = plt.colorbar(image)
    cax.set_label("[km/s]", fontsize=24, rotation=270, labelpad=25)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title(comp_name, loc="right")
    ax.set_title(line_name, loc="left")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 5 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
    if savefig:
        plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/{param}.png", dpi=600)
        plt.close()
    else:
        plt.show()

def render_vel_disp_plot(data, line_center, line_name, wcs, param="G2SIGMA", savefig=False):
    """
    Renders a spaxel map illustrating the number of components fit to 
    each spaxel
    """
    
    # Initialize base aarray
    base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))
    all_sigma = []
    for idx, _ in enumerate(data["XPIX"]):
        vel_disp = ((const.c * (np.abs(data[param][idx]))/line_center).to(u.kilometer / u.second)).value
        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = vel_disp
        all_sigma.append(data["G1SIGMA"][idx])
        all_sigma.append(data["G2SIGMA"][idx])
        all_sigma.append(data["G3SIGMA"][idx])
    
    if "1" in param:
        comp_name = r"$\sigma_{1}$"
    if "2" in param:
        comp_name = r"$\sigma_{2}$"
    if "3" in param:
        comp_name = r"$\sigma_{3}$"
    
    min_disp_1 = np.nanpercentile(((const.c * (np.abs(data["G1SIGMA"]))/line_center).to(u.kilometer / u.second)).value, 1)
    min_disp_2 = np.nanpercentile(((const.c * (np.abs(data["G2SIGMA"]))/line_center).to(u.kilometer / u.second)).value, 1)
    min_disp_3 = np.nanpercentile(((const.c * (np.abs(data["G3SIGMA"]))/line_center).to(u.kilometer / u.second)).value, 1)
    max_disp_1 = np.nanpercentile(((const.c * (np.abs(data["G1SIGMA"]))/line_center).to(u.kilometer / u.second)).value, 99)
    max_disp_2 = np.nanpercentile(((const.c * (np.abs(data["G2SIGMA"]))/line_center).to(u.kilometer / u.second)).value, 99)
    max_disp_3 = np.nanpercentile(((const.c * (np.abs(data["G3SIGMA"]))/line_center).to(u.kilometer / u.second)).value, 99)

    min_disp = np.nanmin([min_disp_1, min_disp_2, min_disp_3])
    max_disp = np.nanmax([max_disp_1, max_disp_2, max_disp_3])
    
    #low_percentile = np.nanpercentile(all_sigma, 1)
    #high_percentile = np.nanpercentile(all_sigma, 99)
    #print(line_center)
    #print(low_percentile, high_percentile)
    
    fig = plt.figure()
    ax = plt.subplot(projection=wcs)
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('gist_heat')
    image = ax.imshow(base_array, vmin=min_disp, vmax=max_disp, cmap=cmap, origin="lower")
    ax.scatter(21, 24, c="white", edgecolors="black", marker="*", s=1000)
    cax = plt.colorbar(image)
    cax.set_label("[km/s]", fontsize=24, rotation=270, labelpad=25)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title(comp_name, loc="right")
    ax.set_title(line_name, loc="left")
    gc_distance = 194.99 * u.Mpc
    scalebar_length = 5 * u.kpc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
    if savefig:
        plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/{param}.png", dpi=600)
        plt.close()
    else:
        plt.show()

def render_vel_disp_rel_vel_scatter(data, line_center, line_name, wcs, savefig=False):
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8)
    prefixes = ["G1", "G2", "G3"]
    rel_vels = []
    disp_vels = []
    amplitudes = []
    for idx, _ in enumerate(data["XPIX"]):
        for prefix in prefixes:
            if not np.isnan(data[prefix + "SIGMA"][idx]):
                if data[prefix + "AMP"][idx] < 1E-4:
                    continue
                rel_vel = ((const.c * (data[prefix + "CEN"][idx] - line_center)/line_center).to(u.kilometer / u.second)).value
                vel_disp = ((const.c * (np.abs(data[prefix + "SIGMA"][idx]))/line_center).to(u.kilometer / u.second)).value
                rel_vels.append(rel_vel)
                disp_vels.append(vel_disp)
                amplitudes.append(data[prefix + "AMP"][idx])
    ax.scatter(rel_vels, disp_vels, c=amplitudes, cmap="viridis", s=5, norm=LogNorm())
    ax.set_xlabel(r"$v$" + " [km/s]")
    ax.set_ylabel(r"$\sigma$" + " [km/s]")
    ax.set_title(line_name, loc="left")
    if savefig:
        plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/vel_disp_rel_vel.png", dpi=600)
        plt.close()
    else:
        plt.show()

line_dict = {"[NeV]": [4, "long", [24.1, 24.5], 24.316],
             "[NeV]_14": [3, "medium", [14.28, 14.36], [14.23, 14.40], 14.3217, [0.03, 0.1, 0.0125, 0.02]],
             "[H_2_S_1]": [3, "long", [16.99, 17.06], [16.98, 17.07], 17.03, [0.03, 0.1, 0.0125, 0.02]], 
             "[H_2_S_3]": [2, "medium", [9.63, 9.71], [9.57, 9.8], 9.6654, [0.03, 0.1, 0.0125, 0.02]], 
             "[OIV]": [4, "long", [25.62, 26.08], 25.88], 
             "MgV": [1, "medium", [5.57, 5.656], 5.609],
             "[NeIII]": [3, "long", [15.5251, 15.6051], [15.38, 15.67], 15.5551, [0.03, 0.1, 0.0125, 0.02]],
             "[SIV]": [2, "long", [10.3, 10.63], [10.00, 10.87], 10.509, [0.03, 0.1, 0.0125, 0.02]],
             "[NeII]": [3, "short", [12.77, 12.84], [12.73, 12.89], 12.813550, [0.03, 0.1, 0.0125, 0.02]],
             "[FeII]": [1, "short", [5.29, 5.37], [5.25, 5.43], 5.3396, [0.01, 0.5, 0.0125, 0.02]]}


line_name = "[NeII]"
data = ascii.read(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/CEN_SNR_fit.dat", format="ipac")  
#data = ascii.read(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/fit.dat", format="ipac")  
spec_obj = CubeSpec("./../", "param_files", "IR23128-S_6_single_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN6", mode="AGN")
south_cubes = spec_obj.perform_single_extraction_custom()
south_data = np.nanmedian(south_cubes[-1].cube_before, axis=0)
wcs = south_cubes[-1].wcs.celestial
plt.style.use('default')
pltparams = PlotParams(palatte="dark", scaling="paper")
keys = list(data.columns)[2:]
for key in keys:
    #print(key)
    if "AMP" in key:
        render_amplitude_plot(data, line_name, wcs, param=key, savefig=True)
    if "CEN" in key:
        render_rel_vel_plot(data, line_dict[line_name][4], line_name, wcs, param=key, savefig=True)
    if "SIGMA" in key:
        render_vel_disp_plot(data, line_dict[line_name][4], line_name, wcs, param=key, savefig=True)

render_vel_disp_rel_vel_scatter(data, line_dict[line_name][4], line_name, wcs, savefig=True)
render_multicomponent_plot(data, wcs, savefig=True)
render_totflux_plot(data, wcs, savefig=True)