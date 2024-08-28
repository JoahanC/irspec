import numpy as np 
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
from astropy.io import ascii
from astropy.visualization import ZScaleInterval


from plotparams import PlotParams


pprams = PlotParams(palatte='dark', scaling="presentation")
#spaxel_arr1 = np.loadtxt("./../diagnostic_plots/spaxel_maps/north_neiii.txt", dtype=float)
#data = ascii.read("./../diagnostic_plots/spaxel_maps/S_[NeIII]_1_1.dat", format="ipac")  
#data = ascii.read("./../diagnostic_plots/spaxel_maps/S_H2_S3_1_1.dat", format="ipac")  

data = ascii.read("./../diagnostic_plots/spaxel_maps/S_[NeV]_14_1_1.dat", format="ipac")  


#base_array = np.zeros((45, 45))
base_array = np.zeros((49, 49))
print(base_array)
for idx, _ in enumerate(data["XPIX"]):
    #if np.abs(data["LINEC"][idx]) > 0.004393:
        rel_vel = ((const.c * (data["LINEC"][idx])/15.555100).to(u.kilometer / u.second)).value
        base_array[data["XPIX"][idx]][data["YPIX"][idx]] = rel_vel
        #base_array[data["XPIX"][idx]][data["YPIX"][idx]] = data["FWHM"][idx] / 100
    #else:
    #    base_array[data["XPIX"][idx]][data["YPIX"][idx]] = np.nan

print(base_array)

#norm=LogNorm(0.01, 0.05), 
z = ZScaleInterval()
z1,z2 = z.get_limits(base_array)
fig, ax = plt.subplots()
image = ax.imshow(base_array, norm=Normalize(-500, 500), origin="lower", cmap="plasma")
ax.set_xlabel("XPIX")
ax.set_ylabel("YPIX")
ax.set_title(r"South: [NeV]$_{14\mu m}$", fontsize=24, loc="right")

#ax.set_title(r"South: [NeV]$_{14\mu m}$", fontsize=24, loc="right")
colorbar = plt.colorbar(image)
colorbar.set_label(r"Relative Veloctiy (km/s)", rotation=270, labelpad=20)

#plt.show()
#colorbar.set_label(r"Flux (Wm$^{-2}$)", fontsize=18, rotation=270, labelpad=25)
plt.savefig("./../diagnostic_plots/spaxel_pres/nev_14_south_relvel.pdf", dpi=600, bbox_inches="tight")