from cubespec import CubeSpec
import numpy as np
import matplotlib.pyplot as plt

datasets = []

entry_idxs = [["S", "6", "AGN"], ["N", "12", "SB"], ["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"]]
    

for entry in entry_idxs:
    spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
    datasets.append(spec_obj.relative_velocities(15.56))


fig, ax = plt.subplots()



for dataset in datasets:
    
    min_dif_array = np.absolute(dataset["relvel"]-2500)
    max_dif_array = np.absolute(dataset["relvel"]+2500)
    low_idx = min_dif_array.argmin()
    high_idx = max_dif_array.argmin() 
    xvals = dataset["relvel"][low_idx:high_idx]
    yvals = dataset["flux"][low_idx:high_idx]
    minval = np.min(yvals)
    newflux = (yvals - minval) / np.max(yvals - minval)
    ax.plot(xvals, newflux)
ax.set_xlim(-2500, 2500)
plt.show()