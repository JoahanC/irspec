import numpy as np 
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
from astropy.io import ascii
from astropy.visualization import ZScaleInterval


from plotparams import PlotParams

core = "S"
name = "[NeV]_14"
npoly = 1
ngauss = 1
param = "LINEC"
pprams = PlotParams(palatte='dark', scaling="presentation")
data = ascii.read(f"./../diagnostic_plots/spaxel_maps/{core}_{name}_{npoly}_{ngauss}.dat", format="ipac")  


#base_array = np.zeros((45, 45))
print(np.max(data["XPIX"]))
base_array = np.zeros((np.max(data["XPIX"]) + 1, np.max(data["YPIX"]) + 1))
for idx, _ in enumerate(data["XPIX"]):
    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = data[param][idx]
    
    """#if np.abs(data["LINEC"][idx]) > 0.004393:
        #rel_vel = ((const.c * (data["LINEC"][idx])/15.555100).to(u.kilometer / u.second)).value
        #base_array[data["XPIX"][idx]][data["YPIX"][idx]] = rel_vel
        print(data["BLUE"][idx])
        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = data["BLUE"][idx]
    #else:
    #    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = np.nan"""

#norm=LogNorm(0.01, 0.05), 
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
#plt.savefig("./../diagnostic_plots/spaxel_pres/h3_south_fwhm.pdf", dpi=1200, bbox_inches="tight")