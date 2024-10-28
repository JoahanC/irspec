import numpy as np 
import os
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
from astropy.io import ascii
from astropy.visualization import ZScaleInterval


from plotparams import PlotParams

core = "S"
name = "[H_2_S_3]"
line_center = 9.6649
type = "triple"
npoly = 1
if type == "double":
    ngauss = 2
if type == "triple":
    ngauss = 3
pprams = PlotParams(palatte='dark', scaling="presentation")
data = ascii.read(f"./../diagnostic_plots/spaxel_maps/{core}_{name}_{npoly}_{ngauss}_multi_{type}.dat", format="ipac")  
base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))

if core not in os.listdir("./../diagnostic_plots/multicomponent_maps"):
    os.mkdir(f"./../diagnostic_plots/multicomponent_maps/{core}")
if name not in os.listdir(f"./../diagnostic_plots/multicomponent_maps/{core}"):
    os.mkdir(f"./../diagnostic_plots/multicomponent_maps/{core}/{name}")
if type not in os.listdir(f"./../diagnostic_plots/multicomponent_maps/{core}/{name}"):
    os.mkdir(f"./../diagnostic_plots/multicomponent_maps/{core}/{name}/{type}")

if type == "double":
    keys = ["FLUXT", "AMP1", "DLINE1", "FWHM1", "AMP2", "DLINE2", "FWHM2"]
if type == "triple":
    keys = ["FLUXT", "AMP1", "DLINE1", "FWHM1", "AMP2", "DLINE2", "FWHM2", "AMP3", "DLINE3", "FWHM3"]
for param in keys:
    
    if "FWHM" in param:
        fwhm_data = data[param]
        fwhm_data[fwhm_data > 1] = np.nan
        for idx, _ in enumerate(data["XPIX"]):
            base_array[data["XPIX"][idx]][data["YPIX"][idx]] = fwhm_data[idx]
    if "LINE" in param:
        for idx, _ in enumerate(data["XPIX"]):
            rel_vel = ((const.c * (data[param][idx] - line_center)/line_center).to(u.kilometer / u.second)).value
            print(data[param][idx])
            base_array[data["XPIX"][idx]][data["YPIX"][idx]] = rel_vel
    if "AMP" in param:
        for idx, _ in enumerate(data["XPIX"]):
            base_array[data["XPIX"][idx]][data["YPIX"][idx]] = data[param][idx]
    
        
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8)
    image = ax.imshow(base_array, origin="lower", cmap="plasma")
    ax.set_xlabel("XPIX")
    ax.set_ylabel("YPIX")
    
    
    if type == "triple":
        if "1" in param:
            ax.set_title(f"{core}: {name} Center", fontsize=24, loc="left")
        if "2" in param:
            ax.set_title(f"{core}: {name} Blue", fontsize=24, loc="left")
        if "3" in param:
            ax.set_title(f"{core}: {name} Red", fontsize=24, loc="left")
    if type == "double":
        if "1" in param:
            ax.set_title(f"{core}: {name} Blue", fontsize=24, loc="left")
        if "2" in param:
            ax.set_title(f"{core}: {name} Red", fontsize=24, loc="left")
    
    if "AMP" in param:
        colorbar = plt.colorbar(image, label="Jy")
        ax.set_title(f"Amplitude", fontsize=24, loc="right")
    if "LINE" in param:
        colorbar = plt.colorbar(image, label="km/s")
        ax.set_title(f"Relative Velocity", fontsize=24, loc="right")
    if "FWHM" in param:
        colorbar = plt.colorbar(image, label="micron")
        ax.set_title(f"FWHM", fontsize=24, loc="right")
    
    plt.savefig(f"./../diagnostic_plots/multicomponent_maps/{core}/{name}/{type}/{param}_spaxel_map.png", dpi=600, bbox_inches="tight")
    plt.close()
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8)
    ax.hist(data[param], bins=100)
    ax.set_title(f"Param: {param}", fontsize=24, loc="right")
    if type == "triple":
        if "1" in param:
            ax.set_title(f"{core}: {name} Center", fontsize=24, loc="left")
        if "2" in param:
            ax.set_title(f"{core}: {name} Blue", fontsize=24, loc="left")
        if "3" in param:
            ax.set_title(f"{core}: {name} Red", fontsize=24, loc="left")
    if type == "double":
        if "1" in param:
            ax.set_title(f"{core}: {name} Blue", fontsize=24, loc="left")
        if "2" in param:
            ax.set_title(f"{core}: {name} Red", fontsize=24, loc="left")
    plt.savefig(f"./../diagnostic_plots/multicomponent_maps/{core}/{name}/{type}/{param}_hist.pdf", dpi=600, bbox_inches="tight")
    plt.close()
    
    """#if np.abs(data["LINEC"][idx]) > 0.004393:
        #rel_vel = ((const.c * (data["LINEC"][idx])/15.555100).to(u.kilometer / u.second)).value
        #base_array[data["XPIX"][idx]][data["YPIX"][idx]] = rel_vel
        print(data["BLUE"][idx])
        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = data["BLUE"][idx]
    #else:
    #    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = np.nan"""

"""#norm=LogNorm(0.01, 0.05), 
#z = ZScaleInterval()
#z1,z2 = z.get_limits(base_array)
fig, ax = plt.subplots()
image = ax.imshow(base_array, origin="lower", cmap="plasma")
ax.set_xlabel("XPIX")
ax.set_ylabel("YPIX")
ax.set_title(r"South: [NeIII] BLUE", fontsize=24, loc="right")

#ax.set_title(r"South: [NeV]$_{14\mu m}$", fontsize=24, loc="right")
colorbar = plt.colorbar(image)
plt.show()
colorbar.set_label(r"Flux (Wm$^{-2}$)", fontsize=18, rotation=270, labelpad=25)
#plt.savefig("./../diagnostic_plots/spaxel_pres/h3_south_fwhm.pdf", dpi=1200, bbox_inches="tight")"""