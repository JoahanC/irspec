import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize, LogNorm
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.io.fits import getdata
import matplotlib.patches as mpatches

from cubespec import CubeSpec


"""z = ZScaleInterval()

filename = "./../archival_data/miri_images/IR23128-N_miri_444/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"
data_array = getdata(filename=filename)
#hdu = fits.open(filename)[0]
#wcs = WCS(hdu.header)
z1, z2 = z.get_limits(data_array)
fig, ax = plt.subplots()
#ax = plt.subplot(projection=wcs)
z1 = 0
z2 = 120
z1 = 5
z2 = 80
#print(hdu.data)
#ax.imshow(data_array, origin='lower')
new_data = data_array[200:350,200:350]
#ax.imshow(data_array, origin="lower", cmap="plasma", norm=LogNorm(z1, z2))
#ax.set_xlim(200, 350)
#ax.set_ylim(200, 350)
z1, z2 = z.get_limits(new_data)
print(z1, z2)
z1=1
z2=150
ax.imshow(new_data, origin="lower", cmap="plasma", norm=LogNorm(z1, z2))
plt.show()"""

filename = "./../creta_output/extractions/IR23128-S/0/IR23128-S_SingleExt_r1.5as.csv"
filename = "./../IR23128-N_SingleExt_r1.5as.csv"

spec_dict = {"w": [], "f": [], "f_unc": [], "DQ": []}
readlines = False
stop_string = "Wave,Band_name,Flux_ap,Err_ap,R_ap,Flux_ap_st,"

err_set = []
with open(filename, 'r') as csvfile:
    for line in csvfile.readlines():
        # Ignore up to this line
        if stop_string in line:
            readlines = True 
            continue
        # Record values
        if readlines:
            vals = line.split(sep=",")
            if vals[0].strip() == "27.906999649014324":
                break
            if vals[5].strip() == "":
                continue
            spec_dict["w"].append(float(vals[0]))
            spec_dict["f"].append(float(vals[5]))
            spec_dict["f_unc"].append(float(vals[6]))
            spec_dict["DQ"].append(float(vals[7]))
            err_set.append(float(vals[7].strip()))
err_set = set(err_set)

colors = []
wavelengths = np.array(spec_dict["w"])
fluxes = np.array(spec_dict["f"])
for idx, flag in enumerate(spec_dict["DQ"]):
    if flag == 0.0:
        colors.append("g")
    if flag == 513.0:
        colors.append("r")
z = 0.044601
fig, ax = plt.subplots()
ax.scatter(wavelengths * (1 - z), fluxes, c=colors, s=2)


spec_obj = CubeSpec("./../", "param_files", "IR23128-N_3_single_param.txt", "input_data/IR23128-N/", redshift=0.044601, fit_dirname="IR23128N_SB3", mode="SB")
af = spec_obj.open_asdf()
wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
flux = np.asarray(af['cafefit']['obsspec']['flux'])
flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
spec_dict = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}

ax.scatter(wave, flux + 0.5, c="black", s=2)
black_patch = mpatches.Patch(color='black', label='CAFE Spectrum')
green_patch = mpatches.Patch(color='green', label='Error Code: 0')
red_patch = mpatches.Patch(color='red', label='Error Code: 513')
ax.set_xlabel("Rest Wavelength (micron)")
ax.set_ylabel(r"$f_{\nu}$ (Relative)")
plt.legend(handles=[black_patch, red_patch, green_patch])

#plt.show()
plt.close()




