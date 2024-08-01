import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits

original_file = "./../creta_output/extractions/IR23128-S/5/IR23128-S_SingleExt_r1.75as.fits"
new_temp_file = "./../original_IR23128-N_SingleExt_r1.8as.fits"

hdul1 = fits.open(original_file)
data_array1 = hdul1[1].data[0]
wave01 = data_array1[0]
flux1 = data_array1[2]
flags1 = data_array1[7]    
hdul1.close()



fig, ax = plt.subplots()

mod_red_wave = []
mod_red_flux = []
mod_green_wave = []
mod_green_flux = []
orig_red_wave = []
orig_red_flux = []
orig_green_wave = []
orig_green_flux = []
for idx, wave in enumerate(wave01):
    if flags1[idx] == 513.0:
        mod_red_wave.append(wave)
        mod_red_flux.append(flux1[idx])
    if flags1[idx] == 0.0:
        mod_green_wave.append(wave)
        mod_green_flux.append(flux1[idx])

"""for idx, wave in enumerate(wave02):
    if flags2[idx] == 513.0:
        orig_red_wave.append(wave)
        orig_red_flux.append(flux2[idx] * 100)
    if flags2[idx] == 0.0:
        orig_green_wave.append(wave)
        orig_green_flux.append(flux2[idx] * 100)"""

ax.scatter(mod_red_wave, mod_red_flux, c="red")
ax.scatter(mod_green_wave, mod_green_flux, c="green")
#ax.scatter(orig_red_wave, orig_red_flux, c="red")
#ax.scatter(orig_green_wave, orig_green_flux, c="green")
ax.set_yscale("log")
plt.show()