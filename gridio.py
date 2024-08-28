import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.wcs import WCS

import astropy.io.fits as fits 
import matplotlib.animation as animation
from cubespec import CubeSpec
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()

from plotparams import PlotParams
from matplotlib.gridspec import GridSpec
pprams = PlotParams(palatte="dark", scaling="presentation")


"""hdul = fits.open("./../creta_output/extractions/IR23128-S/11/IR23128-S_GridExt_5x5_s1.0as.fits")
print(hdul.info())
keys = []
data = hdul[1].data
for i in range(8):
    keys.append(hdul[1].header[f"TTYPE{i + 1}"])

fig, ax = plt.subplots()
for region in data: 
    ax.plot(region[0], region[2])
ax.set_yscale("log")
plt.show()"""

"""hdul = fits.open("./../input_data/IR23128-S_grid/Level3_ch4-long_s3d.fits")
print(hdul.info())
print(repr(hdul[1].header))
wave0 = hdul[1].header["CRVAL3"]
dwave = hdul[1].header["CDELT3"]
steps = np.arange(717)
deltawave = steps * dwave
wavelengths = wave0 + deltawave

data = hdul[1].data
print(np.shape(data))



fig, ax = plt.subplots()

for i in range(1):
    for j in range(1):
        spectrum_1 = data[:, 10, 10]
        ax.plot(wavelengths, spectrum_1)
        print(i ,j)
ax.set_yscale("log")
plt.show()"""

"""hdul = fits.open("./../input_data/IR23128-N_grid/iras23128-n_ch1-medium_s3d.fits")
print(hdul.info())
print(repr(hdul[1].header))
wave0 = hdul[1].header["CRVAL3"]
dwave = hdul[1].header["CDELT3"]
steps = np.arange(1213)
deltawave = steps * dwave
wavelengths = wave0 + deltawave

data = hdul[1].data
print(np.shape(data))

print(wavelengths[550], wavelengths[800])


fig, ax = plt.subplots()

for i in range(49):
    for j in range(47):
        spectrum_1 = data[550:800, i, j]
        ax.plot(wavelengths[550:800], spectrum_1)
        print(i ,j)
ax.set_yscale("log")
plt.show()"""

"""hdul = fits.open("./../input_data/IR23128-N_grid/iras23128-n_ch2-medium_s3d.fits")
print(hdul.info())
print(repr(hdul[1].header))
wave0 = hdul[1].header["CRVAL3"]
dwave = hdul[1].header["CDELT3"]
steps = np.arange(1124)
deltawave = steps * dwave
wavelengths = wave0 + deltawave

data = hdul[1].data
print(np.shape(data))



print(wavelengths[755], wavelengths[770])


fig, ax = plt.subplots()

for i in range(45):
    for j in range(45):
        spectrum_1 = data[755:770, i, j]
        ax.plot(wavelengths[755:770], spectrum_1)
        print(i ,j)
ax.set_yscale("log")
plt.show()

#print(spectrum_1)"""


"""fig = plt.figure()
fig.suptitle("Controlling spacing around and between subplots")
#fig, ax = plt.subplots()

gs1 = GridSpec(2, 2, left=0.05, right=0.95, wspace=0.05)
ax1 = fig.add_subplot(gs1[0,0])
ax2 = fig.add_subplot(gs1[0,1])
ax = fig.add_subplot(gs1[1,:])

spec_obj = CubeSpec("./../", "param_files", "IR23128-S_0_single_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN0", mode="AGN")

dataset = spec_obj.recall_data()

t = dataset["wave"]
z = dataset["flux"]
vl = ax.axvline(t[0], ls='-', color='yellow', lw=1)

scat = ax.scatter(t[0], z[0], c="b", s=5, label=f'void')
ax.set(xlim=[np.min(t), np.max(t)], ylim=[np.min(z), np.max(z)], xlabel='Time [s]', ylabel='Z [m]')
ax.set_yscale("log")
ax.legend()
print(len(t))


def update(frame):
    # for each frame, update the data stored on each artist.
    x = t[:frame * 5]
    y = z[:frame * 5]
    if frame != 0:
        vl.set_xdata([x[-1],x[-1]])
    # update the scatter plot:
    data = np.stack([x, y]).T
    scat.set_offsets(data)

    return scat


ani = animation.FuncAnimation(fig=fig, func=update, frames=len(t)//5, interval=1)
#ani.save('orbita.gif', writer='imagemagick', fps=30)
plt.show()"""

redshift = 0.044601
datacube_s = f"./../input_data/IR23128-S/Level3_ch2-short_s3d.fits"
datacube_n = f"./../input_data/IR23128-N/iras23128-n_ch2-short_s3d.fits"

hdul_n = fits.open(datacube_n)
data_n = hdul_n[1].data
wcs_n = WCS(hdul_n[1].header)
wave0_n = hdul_n[1].header["CRVAL3"]
dwave_n = hdul_n[1].header["CDELT3"]

steps_n = np.arange(len(data_n))
deltawave_n = steps_n * dwave_n
wavelengths_n = (wave0_n + deltawave_n) / (1 + redshift)
fluxes_n = data_n[:, 26, 21]

hdul_s = fits.open(datacube_s)
data_s = hdul_s[1].data
wave0_s = hdul_s[1].header["CRVAL3"]
dwave_s = hdul_s[1].header["CDELT3"]

steps_s = np.arange(len(data_s))
deltawave_s = steps_s * dwave_s
wavelengths_s = (wave0_s + deltawave_s) / (1 + redshift)
fluxes_s = data_s[:, 25, 20]

fig = plt.figure()
fig.set_size_inches(16, 10)
#fig, ax = plt.subplots()

t_n = np.copy(wavelengths_n)
z_n = np.copy(fluxes_n)
cubes_n = np.copy(data_n)

t_s = np.copy(wavelengths_s)
z_s = np.copy(fluxes_s)
cubes_s = np.copy(data_s)

gs1 = GridSpec(2, 2, left=0.1, right=0.9, wspace=0.25, hspace=0.15)
ax1 = fig.add_subplot(gs1[0,0])
ax2 = fig.add_subplot(gs1[0,1])

#ax1 = plt.subplot(gs1[0,0], projection=wcs_n, slices=('x', 'y', 0))
#ax2 = plt.subplot(gs1[0,1], projection=wcs_s, slices=('x', 'y', 0))

ax = fig.add_subplot(gs1[1,:])

z1_n,z2_n = zscale.get_limits(cubes_n[0])
z1_s,z2_s = zscale.get_limits(cubes_s[0])

scatterp_n = ax.scatter(t_n[0], z_n[0], c="lime", s=5)
scatterp_s = ax.scatter(t_s[0], z_s[0], c="white", s=5)

cubeplot = ax1.imshow(cubes_n[0], cmap="plasma", origin="lower")
cubeplot2 = ax2.imshow(cubes_s[0], cmap="plasma", origin="lower")

title = ax.text(np.average(t_n)*0.97, np.max(z_s)*90, f"Displayed Wavelength: {round(wavelengths_n[0], 3)} microns", fontdict={'size': 20})

plots = [scatterp_n, scatterp_s, cubeplot, cubeplot2, title]
vl = ax.axvline(t_n[0], ls='-', color='red', lw=1)

ax1.set(xlabel='XPIX', ylabel='YPIX')
ax2.set(xlabel='XPIX', ylabel='YPIX')
ax1.set_title("Northern Core", loc="left", fontsize=24)
ax2.set_title("Southern Core", loc="left", fontsize=24)
min_y = np.min([np.min(z_n), np.min(z_s)])
max_y = np.max([np.max(z_n), np.max(z_s)])
ax.set(xlim=[np.min(t_n), np.max(t_n)], ylim=[min_y*0.9, max_y*1.2], xlabel=r'$\lambda_{rest}$ (micron)', ylabel=r'$f_{\nu}$ (Jy)')
ax.set_yscale("log")
ax.legend()

circle1 = plt.Circle((21, 26), 3, color='lime', lw=3, fill=False)
circle2 = plt.Circle((20, 25), 3, color='white', lw=3, fill=False)

ax1.add_patch(circle1)
ax2.add_patch(circle2)



print(fluxes_n[0:5])
print(fluxes_s[0:5])

def update(frame):
    # for each frame, update the data stored on each artist.
    x_n = t_n[:frame * 5]
    y_n = z_n[:frame * 5]
    
    x_s = t_s[:frame * 5]
    y_s = z_s[:frame * 5]
    
    cubes2_n = cubes_n[frame*5]
    cubes2_s = cubes_s[frame*5]
    
    # update the scatter plot:
    data1_n = np.stack([x_n, y_n]).T
    data1_s = np.stack([x_s, y_s]).T
    
    
    
    if frame != 0:
        vl.set_xdata([x_n[-1],x_n[-1]])
        plots[4].set_text(f"Displayed Wavelength: {round(x_n[-1], 3)} microns")
    
    
    
    plots[0].set_offsets(data1_n)
    plots[1].set_offsets(data1_s)
    plots[2].set_data(cubes2_n)
    plots[3].set_data(cubes2_s)
    
    
    print(frame, plots)
    #return scatterp
    return plots

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(wavelengths_n)//5, interval=30)
#ani.save('orbita.gif', writer='imagemagick', fps=30)
plt.show()
