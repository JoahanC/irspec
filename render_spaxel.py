import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm

from plotparams import PlotParams

pprams = PlotParams(palatte='dark', scaling="presentation")
spaxel_arr1 = np.loadtxt("./../diagnostic_plots/spaxel_maps/north_neiii.txt", dtype=float)
print(spaxel_arr1[20][20])
fig, ax = plt.subplots()
image = ax.imshow(spaxel_arr1, norm=LogNorm(), origin="lower", cmap="plasma")
ax.set_xlabel("XPIX")
ax.set_ylabel("YPIX")
ax.set_title(r"North: [NeIII]", fontsize=24, loc="right")
colorbar = plt.colorbar(image)
colorbar.set_label(r"Flux (Wm$^{-2}$)", fontsize=18, rotation=270, labelpad=25)
plt.savefig("./../diagnostic_plots/spaxel_pres/ne3_north.pdf", dpi=1200, bbox_inches="tight")