import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import astropy.io.fits as fits 
from astropy.wcs import WCS
from astropy.table import Table


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
blue_linec = []
core_linec = []
red_linec = []


waverange = [15.48, 15.53]

for i in range(np.shape(data)[1]):
    for j in range(np.shape(data)[2]):
        x_pix.append(i)
        y_pix.append(j) 
        if data_quality[i, j] == 513:
            blue_linec.append(np.nan)
            #core_linec.append(np.nan)
            #red_linec.append(np.nan)
            continue
        full_spec_1 = data[:, i, j]
        
        
        
        loc = LocalFit(wavelengths, full_spec_1, waverange, wave_c, name)
        if i == 13 and j == 32:
            loc.render_data()
        #loc.neon3_fit()
        
        
        if loc.neon3_fit() == None:
            blue_linec.append(np.nan)
            core_linec.append(np.nan)
            red_linec.append(np.nan)
        else:
            print(loc.popt)
            blue_linec.append(loc.popt[2])
            core_linec.append(loc.popt[5])
            red_linec.append(loc.popt[8])

        print(f"SPAXEL: {i}, {j}")
        


fitparams = Table([x_pix, y_pix, blue_linec, core_linec, red_linec], names=("XPIX", "YPIX", "BLUE", "CORE", "RED"))#
fitparams.write(f"./../diagnostic_plots/spaxel_maps/neiii_fluxes.dat", format="ipac", overwrite=True)

"""np.savetxt(f"./../diagnostic_plots/spaxel_maps/{core}_{name}.txt", fluxes)

pprams = PlotParams(palatte="dark", scaling="presentation")

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