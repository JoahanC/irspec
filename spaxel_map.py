import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import astropy.io.fits as fits 
from astropy.wcs import WCS
from astropy.table import Table
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval
z = ZScaleInterval()

from plotparams import PlotParams


from localfit import LocalFit


core="S"
name = "[NeV]_14"

# NeV
#channel = 4
#subchannel = "long"
#waverange = [24.1, 24.5]
#wave_c = 24.316

# NeV 14micron
#channel = 3
#subchannel = "medium"
#waverange = [14.28, 14.35]
#wave_c = 14.319
#true_wave_c = 14.3217

# H2 S3
#channel = 2
#subchannel = "medium"
#waverange = [9.64, 9.69]
#wave_c = 9.667
#true_wave_c = 9.6649

# OIV
#channel = 4
#subchannel = "long"
#waverange = [25.62, 26.08]
#wave_c = 25.88

# MgV FIX
#channel = 1
#subchannel = "medium"
#waverange = [5.57, 5.656]
#wave_c = 5.609

# NeIII
channel = 3
subchannel = "long"
waverange = [15.38, 15.67]
wave_c = 15.555100

# NeII
#channel = 3
#subchannel = "medium"
#waverange = [12.77, 12.9]
#wave_c = 12.81
#true_wave_c = 

redshift = 0.044601
if core == "N":
    datacube = f"./../input_data/IR23128-N_grid/iras23128-n_ch{channel}-{subchannel}_s3d.fits"
if core == "S":
    datacube = f"./../input_data/IR23128-S/Level3_ch{channel}-{subchannel}_s3d.fits"


hdul = fits.open(datacube)
data = hdul[1].data
wave0 = hdul[1].header["CRVAL3"]
dwave = hdul[1].header["CDELT3"]

steps = np.arange(len(data))
deltawave = steps * dwave
wavelengths = (wave0 + deltawave) / (1 + redshift)

data_quality = hdul[3].data[0]
fluxes = np.zeros((np.shape(data)[1], np.shape(data)[2]))
#print(np.shape(data))

print(wavelengths[0], wavelengths[1])

x_pix = []
y_pix = []
fluxes_t = []
line_center = []
fwhm = []

"""
pixels = [(13, 24), (15, 28), (16,31), (6,27)]
blue_pixels = [(28, 15),  (27,6), (33, 12)]
red_pixels = [(27, 40), (35, 29), (25,32), (23,38), (33, 12)]
red_pixels = [(21, 26),  (30,30), (29, 35)]


pprams = PlotParams(scaling="presentation")
for pixel in blue_pixels:
    full_spec_1 = data[:, pixel[0], pixel[1]]
    full_spec_2 = data[:, pixel[0]+1, pixel[1]]
    full_spec_3 = data[:, pixel[0], pixel[1]+1]
    full_spec_4 = data[:, pixel[0]+1, pixel[1]+1]
    full_spec = full_spec_1 + full_spec_2 + full_spec_3 + full_spec_4
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.plot(wavelengths, full_spec / np.max(full_spec), color="white", lw=5)
    ax.set_xlim(15.475, 15.6)
    #ax.set_yscale("log")
    ax.axvline(15.553, color="yellow", alpha=1)
    #ax.set_ylabel("Relative In")
    plt.xticks([])
    plt.yticks([])
    ax.grid()
    plt.savefig(f"./line_profile_blue_{pixel[0]}_{pixel[1]}.pdf", dpi=600, bbox_inches="tight")
    plt.close()

for idx, pixel in enumerate(red_pixels):
    full_spec_1 = data[:, pixel[0], pixel[1]]
    full_spec_2 = data[:, pixel[0]+1, pixel[1]]
    full_spec_3 = data[:, pixel[0], pixel[1]+1]
    full_spec_4 = data[:, pixel[0]+1, pixel[1]+1]
    full_spec = full_spec_1 + full_spec_2 + full_spec_3 + full_spec_4
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.plot(wavelengths, full_spec / np.max(full_spec), color="white", lw=5)
    ax.set_xlim(15.475, 15.6)
    #ax.set_yscale("log")
    ax.axvline(15.553, color="yellow", alpha=1)
    plt.xticks([])
    plt.yticks([])
    #ax.set_ylabel("Relative In")
    ax.grid()
    plt.savefig(f"./line_profile_red_{pixel[0]}_{pixel[1]}.pdf", dpi=600, bbox_inches="tight")
    plt.close()

blue_dif_array = np.absolute(wavelengths - (15.533))
difference = 15.5526 - 15.533
red_dif_array = np.absolute(wavelengths - (15.5526 + difference))
blue_idx = blue_dif_array.argmin()
red_idx = red_dif_array.argmin()

print(blue_idx, red_idx)
blue_cube = data[blue_idx]
red_cube = data[red_idx]

z1, z2 = z.get_limits(blue_cube)
fig, ax = plt.subplots()
ax.imshow(blue_cube, norm=LogNorm(1, z2), origin="lower", cmap="plasma")
#ax.set_title(r"15.535 $\mu m$", fontsize=24)
ax.set_title(f"-{round(difference, 4)}" + r" $\mu m$", fontsize=24)
for pixel in blue_pixels:
    circle=plt.Circle((pixel[1],pixel[0]),2, color="black", fill=False, lw=3)
    ax.add_patch(circle)
ax.set_xlabel("XPIX")
ax.set_ylabel("YPIX")
plt.savefig("blue_cube.pdf", dpi=600, bbox_inches="tight")

z1, z2 = z.get_limits(red_cube)
fig, ax = plt.subplots()

ax.imshow(red_cube, norm=LogNorm(0.5, z2), origin="lower", cmap="plasma")
for pixel in red_pixels:
    circle=plt.Circle((pixel[1],pixel[0]),2, color="black", fill=False, lw=3)
    ax.add_patch(circle)
ax.set_xlabel("XPIX")
ax.set_ylabel("YPIX")
ax.set_title(f"+{round(difference, 4)}" + r" $\mu m$", fontsize=24)
plt.savefig("red_cube.pdf", dpi=600, bbox_inches="tight")"""


"""for i in range(np.shape(data)[1]):
    for j in range(np.shape(data)[2]):
        x_pix.append(i)
        y_pix.append(j) 
        if data_quality[i, j] == 513:
            fluxes_t.append(np.nan)
            line_center.append(np.nan)
            fwhm.append(np.nan)
            continue
        full_spec_1 = data[:, i, j]
        
        
        
        #fig, ax = plt.subplots()
        #ax.plot(wavelengths, full_spec_1)
        #ax.set_yscale("log")
        #plt.show()
        loc = LocalFit(wavelengths, full_spec_1, waverange, wave_c, name)
        #loc.main_fit(npoly=1, ngauss=1, spaxel_fit=True)
        loc.neon3_fit()
        #if i == 15 and j == 30:
        #    loc.render_fit(path="./test.pdf")
        print(loc.popt)

        print(f"SPAXEL: {i}, {j}")
        #print(loc.line_strength.value)
        if loc.line_strength == 0:
            fluxes_t.append(np.nan)
            line_center.append(np.nan)
            fwhm.append(np.nan)
            fluxes[i, j] = 0
        else:
            fluxes_t.append(loc.line_strength.value)
            line_center.append(loc.popt[3] - true_wave_c)
            fwhm.append(loc.popt[3] * 2.355)
            fluxes[i, j] = loc.line_strength.value"""

"""fitparams = Table([x_pix, y_pix, fluxes_t, line_center, fwhm], names=("XPIX", "YPIX", "FLUX", "LINEC", "FWHM"))
fitparams.write(f"./../diagnostic_plots/spaxel_maps/{core}_{name}_{loc.npoly}_{loc.ngauss}.dat", format="ipac", overwrite=True)

np.savetxt(f"./../diagnostic_plots/spaxel_maps/{core}_{name}.txt", fluxes)"""

"""pprams = PlotParams(palatte="dark", scaling="presentation")

fig = plt.figure()
wcs = WCS(hdul[1].header)
ax = plt.subplot(projection=wcs, slices=('x', 'y', 20))
fig.set_size_inches(10, 10)
cbar = ax.imshow(fluxes, cmap="plasma", origin="lower", norm=mpl.colors.LogNorm())
ax.set_xlabel("RA (J2000)")
ax.set_ylabel("Dec (J2000)")
ax.set_title(f"{name} Flux Map")
plt.colorbar(cbar)
plt.savefig(f"./../diagnostic_plots/spaxel_maps/{core}_{name}.pdf", dpi=1200, bbox_inches="tight")"""