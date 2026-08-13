import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from matplotlib.colors import LogNorm, Normalize
from tqdm import tqdm
import astropy.units as u 
import astropy.constants as const
import matplotlib.font_manager as fm
from astropy.io import ascii
from irspec.emission_io import read_line_params, read_line_params2
from irspec import paths


def load_fit(filepath):
    """Loads a saved instance of Spaxelcube"""
    return ascii.read(filepath, format="ipac")

def render_totflux_plot(fitparams, line_dict, name, savefig=False):
    """
    Renders a spaxel map illustrating the number of components fit to 
    each spaxel
    """
    
    # Initialize base aarray
    base_array = np.zeros((np.max(fitparams["XPIX"]) + 1, np.max(fitparams["YPIX"]) + 1))
    
    for idx, _ in enumerate(fitparams["XPIX"]):
        flux_val = 0
        if not np.isnan(fitparams["G1AMP"][idx]):
            
            g1_gamma = np.abs(fitparams["G1SIGMA"][idx]) * 2.355 / fitparams["G1CEN"][idx]
            g1_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (fitparams["G1AMP"][idx] * u.Jy * g1_gamma / (fitparams["G1CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
            flux_val += g1_flux.value
            
            if not np.isnan(fitparams["G2AMP"][idx]):
                g2_gamma = np.abs(fitparams["G2SIGMA"][idx]) * 2.355 / fitparams["G2CEN"][idx]
                g2_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (fitparams["G2AMP"][idx] * u.Jy * g2_gamma / (fitparams["G2CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                flux_val += g2_flux.value
                if not np.isnan(fitparams["G3AMP"][idx]):
                    g3_gamma = np.abs(fitparams["G3SIGMA"][idx]) * 2.355 / fitparams["G3CEN"][idx]
                    g3_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (fitparams["G3AMP"][idx] * u.Jy * g3_gamma / (fitparams["G3CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                    flux_val += g3_flux.value
                    if flux_val < 0:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = 0
                    else:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
                else:

                    if flux_val < 0:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = 0
                    else:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
            else:

                if flux_val < 0:
                    base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = 0
                else:
                    base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
        else:
            base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val

    fig = plt.figure()
    ax = plt.subplot()
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    image = ax.imshow(base_array, cmap=cmap, norm=LogNorm(), origin="lower")
    if line_dict[name][0] == 1:
        ax.scatter(24, 28, c="white", edgecolors="black", marker="*", s=1000)
    if line_dict[name][0] == 2:
        ax.scatter(21, 26, c="white", edgecolors="black", marker="*", s=1000)
    if line_dict[name][0] == 3:
        ax.scatter(21, 25, c="white", edgecolors="black", marker="*", s=1000)
    if line_dict[name][0] == 4:
        ax.scatter(18, 21, c="white", edgecolors="black", marker="*", s=1000)
    cax = plt.colorbar(image)
    cax.set_label(r"[Jy]", fontsize=24, rotation=270, labelpad=25)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title("Total Flux", loc="right")
    title_name = name
    if "H_2_S_1" in name:
        title_name = r"H$_{2}$ S(1)"
    if "H_2_S_3" in name:
        title_name = r"H$_{2}$ S(3)"
    if "14" in name:
        title_name = "[NeV]"+r"$_{14}$"
    ax.set_title(title_name, loc="left")
    #c_distance = 194.99 * u.Mpc
    #calebar_length = 5 * u.kpc
    #scalebar_angle = (scalebar_length / gc_distance).to(
    #    u.deg, equivalencies=u.dimensionless_angles()
    #)
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    #add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
    if savefig:
        outdir = paths.plots_dir() / "dynamic_multicomponent" / name
        outdir.mkdir(parents=True, exist_ok=True)
        plt.savefig(outdir / "flux_white.png", dpi=600)
        plt.close()
    else:
        plt.show()

def calc_totflux(fitparams):
    """
    Calculate the total flux array.
    """
    
    # Initialize base aarray
    base_array = np.zeros((np.max(fitparams["XPIX"]) + 1, np.max(fitparams["YPIX"]) + 1))
    
    for idx, _ in enumerate(fitparams["XPIX"]):
        flux_val = 0
        if not np.isnan(fitparams["G1AMP"][idx]):
            
            g1_gamma = np.abs(fitparams["G1SIGMA"][idx]) * 2.355 / fitparams["G1CEN"][idx]
            g1_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (fitparams["G1AMP"][idx] * u.Jy * g1_gamma / (fitparams["G1CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
            flux_val += g1_flux.value
            
            if not np.isnan(fitparams["G2AMP"][idx]):
                g2_gamma = np.abs(fitparams["G2SIGMA"][idx]) * 2.355 / fitparams["G2CEN"][idx]
                g2_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (fitparams["G2AMP"][idx] * u.Jy * g2_gamma / (fitparams["G2CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                flux_val += g2_flux.value
                if not np.isnan(fitparams["G3AMP"][idx]):
                    g3_gamma = np.abs(fitparams["G3SIGMA"][idx]) * 2.355 / fitparams["G3CEN"][idx]
                    g3_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (fitparams["G3AMP"][idx] * u.Jy * g3_gamma / (fitparams["G3CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                    flux_val += g3_flux.value
                    if flux_val < 0:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = 0
                    else:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
                else:

                    if flux_val < 0:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = 0
                    else:
                        base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
            else:

                if flux_val < 0:
                    base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = 0
                else:
                    base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
        else:
            base_array[fitparams["XPIX"][idx]][fitparams["YPIX"][idx]] = flux_val
    return base_array


def plot_array(base_array, title, log=False):
    fig = plt.figure()
    ax = plt.subplot()
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    if log:
        image = ax.imshow(np.log10(base_array), origin="lower")
    if not log:
        image = ax.imshow(base_array, vmin=0.1, vmax=10, origin="lower")
    cax = plt.colorbar(image)
    cax.set_label(r"[Jy]", fontsize=24, rotation=270, labelpad=25)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title(title, loc="right")
    fontprops = fm.FontProperties(size=24, family='Helvetica')
    plt.show()

def plot_fe2ar2(base_array, title):    
    fig = plt.figure()
    ax = plt.subplot()
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    image = ax.imshow(base_array, vmin=0.1, vmax=2, origin="lower")
    cax = plt.colorbar(image)
    cax.set_label("log([FeII]/[ArII])", fontsize=24, rotation=270, labelpad=25)
    cax.ax.tick_params(labelsize=18)
    ax.set_xlabel("XPIX", fontsize=20)
    ax.set_ylabel("YPIX", fontsize=20)
    #ax.set_title("", loc="right")
    ax.tick_params(axis='both', which='major', labelsize=18)   
    #fontprops = fm.FontProperties(size=24, family='Helvetica')
    plt.savefig(paths.spaxs_dir() / "[FeIIArII]" / "fe2ar2_ratio.png", dpi=400)

def plot_ne3ne2(base_array):
    fig = plt.figure()
    ax = plt.subplot()
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    image = ax.imshow(np.log10(base_array), origin="lower")
    cax = plt.colorbar(image)
    cax.set_label("log([NeIII]/[NeII])", fontsize=24, rotation=270, labelpad=25)
    cax.ax.tick_params(labelsize=18)
    ax.set_xlabel("XPIX", fontsize=20)
    ax.set_ylabel("YPIX", fontsize=20)
    #ax.set_title("[NeIII]/[NeII]", loc="right")
    #fontprops = fm.FontProperties(size=24, family='Helvetica')
    ax.tick_params(axis='both', which='major', labelsize=18)   
    plt.savefig(paths.spaxs_dir() / "[NeIIINeII]" / "ne3ne2_ratio.png", dpi=400)

def plot_ne5ne2(base_array):
    fig = plt.figure()
    ax = plt.subplot()
    fig.set_size_inches(10, 8)
    cmap = plt.get_cmap('plasma')
    image = ax.imshow(np.log10(base_array), vmin=-2.5, vmax=-0.5, origin="lower")
    cax = plt.colorbar(image)
    cax.set_label(r"log([NeV]$_{14}$/[NeII])", fontsize=24, rotation=270, labelpad=25)
    cax.ax.tick_params(labelsize=18)
    ax.set_xlabel("XPIX", fontsize=20)
    ax.set_ylabel("YPIX", fontsize=20)
    #fontprops = fm.FontProperties(size=24, family='Helvetica')
    ax.tick_params(axis='both', which='major', labelsize=18)   
    plt.savefig(paths.spaxs_dir() / "[NeVNeII]" / "ne5ne2_ratio.png", dpi=400)

def neon_scatter(ne5ne2_arr, ne3ne2_arr):
    import matplotlib.patches as patches

    ne5ne2_vals = []
    ne3ne2_vals = []
    for i in range(np.shape(ne5ne2_arr)[1]):
        for j in range(np.shape(ne5ne2_arr)[0]):
            if not np.isnan(ne5ne2_arr[i][j]) and not np.isnan(ne3ne2_arr[i][j]):
                ne3ne2_vals.append(ne3ne2_arr[i][j])
                ne5ne2_vals.append(ne5ne2_arr[i][j])
    fig = plt.figure()
    ax = plt.subplot()
    fig.set_size_inches(10, 8)
    ax.scatter(np.log10(ne3ne2_vals), np.log10(ne5ne2_vals), color="black", s=20, label="Spaxels")
    sfagn_x = [-0.45, -0.15, 0.15, 0.45, 0.75]
    sfagn_y = [-1.1, -1.0, -0.9, -0.75, -0.55]
    ax.scatter(sfagn_x, sfagn_y, color="purple")
    ax.plot(sfagn_x, sfagn_y, color="purple", ls="dashed", label="AGN + SFR")
    rect = patches.Rectangle((-0.5, -1.3), 0.8, 1.1, linewidth=1, edgecolor='green', facecolor='green', alpha=0.5, zorder=1, label="AGN + Shocks")
    ax.add_patch(rect)
    ax.axvspan(xmin=-0.20, xmax=0.9, color="red", alpha=0.3, zorder=0, label="AGN")
    ax.set_xlabel("log([NeIII]/[NeII])", fontsize=20)
    ax.set_xlim(-1.6, 0.8)
    ax.set_ylabel(r"log([NeV]$_{14}$/[NeII])", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)    
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.legend(prop={'size': 18})
    plt.savefig(paths.spaxs_dir() / "[NeVNeII]" / "neon_scatter.png", dpi=400)


def main():
    spaxs = paths.spaxs_dir()
    ne3_path = spaxs / "test" / "[NeIII]" / "twogaussian_raw.dat"
    ne5_path = spaxs / "test" / "[NeV]_14" / "twogaussian_raw.dat"
    ne2_path = spaxs / "test" / "[NeII]" / "twogaussian_raw.dat"
    ne2_fitparams = load_fit(ne2_path)
    ne3_fitparams = load_fit(ne3_path)
    ne5_fitparams = load_fit(ne5_path)

    ne2_arr = calc_totflux(ne2_fitparams)
    ne3_arr = calc_totflux(ne3_fitparams)
    ne5_arr = calc_totflux(ne5_fitparams)
    ne3ne2_arr = ne3_arr / ne2_arr
    ne5ne2_arr = ne5_arr / ne2_arr

    fe2_path = spaxs / "test" / "[FeII]" / "twogaussian_raw.dat"
    ar2_path = spaxs / "test" / "[ArII]" / "twogaussian_raw.dat"
    fe2_fitparams = load_fit(fe2_path)
    ar2_fitparams = load_fit(ar2_path)

    fe2_arr = calc_totflux(fe2_fitparams)
    ar2_arr = calc_totflux(ar2_fitparams)
    fe2ar2_arr = fe2_arr / ar2_arr
    plot_fe2ar2(fe2ar2_arr, "[FeII]/[ArII]")
    plot_ne3ne2(ne3ne2_arr)
    plot_ne5ne2(ne5ne2_arr)


if __name__ == "__main__":
    main()
