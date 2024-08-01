import numpy as np 
from astropy.io import fits
import subprocess

original_file = "./../creta_output/extractions/IR23128-S/6/IR23128-S_SingleExt_r1.7as.fits"
new_temp_file = "./../creta_output/extractions/IR23128-S/1/original_IR23128-S_SingleExt_r1.7as.fits"

subprocess.run(["cp", original_file, new_temp_file])


hdul = fits.open(new_temp_file)
data_array = hdul[1].data[0]


new_flags = np.copy(data_array[7])
for idx, val in enumerate(data_array[7]):
    if val == 513.:
        if not np.isnan(data_array[2][idx]):
            print(data_array[2][idx])
            new_flags[idx] = 0.0
hdul[1].data[0][7] = new_flags
hdul.writeto(original_file, overwrite=True)
hdul.close()