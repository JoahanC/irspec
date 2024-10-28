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


core = "S"
name = "[H_2_S_3]"
type = "triple"
multicomponent = True

line_dict = {"[NeV]": [4, "long", [24.1, 24.5], 24.316],
             "[NeV]_14": [3, "medium", [14.28, 14.35], 14.3217],
             "[H_2_S_3]": [2, "medium", [9.64, 9.69], 9.6649], 
             "[OIV]": [4, "long", [25.62, 26.08], 25.88], 
             "MgV": [1, "medium", [5.57, 5.656], 5.609],
             "[NeIII]": [3, "long", [15.38, 15.67], 15.5551],
             "[NeII]": [3, "medium", [12.77, 12.9], 12.81],
             "[FeII]": [1, "short", [5.30, 5.41], 5.3396]}

redshift = 0.044601
if core == "N":
    datacube = f"./../input_data/IR23128-N/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"
if core == "S":
    datacube = f"./../input_data/IR23128-S/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"


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

"""fig, ax = plt.subplots()
ax.imshow(data_quality)
plt.show()"""

print(wavelengths[0], wavelengths[1])

x_pix = []
y_pix = []
fluxes_t = []
line_center = []
fwhm = []
fluxes_amp_1 = []
line_center_1 = []
fwhm_1 = []
fluxes_amp_2 = []
line_center_2 = []
fwhm_2 = []
fluxes_amp_3 = []
line_center_3 = []
fwhm_3 = []

for i in range(np.shape(data)[1]):
    for j in range(np.shape(data)[2]):
        
        x_pix.append(i)
        y_pix.append(j) 
        if data_quality[i, j] == 513:
            if multicomponent:
                fluxes_t.append(np.nan)
                fluxes_amp_1.append(np.nan)
                line_center_1.append(np.nan)
                fwhm_1.append(np.nan)
                fluxes_amp_2.append(np.nan)
                line_center_2.append(np.nan)
                fwhm_2.append(np.nan)
                if type == "triple":
                    fluxes_amp_3.append(np.nan)
                    line_center_3.append(np.nan)
                    fwhm_3.append(np.nan)
                continue
            else:
                fluxes_t.append(np.nan)
                line_center.append(np.nan)
                fwhm.append(np.nan)
                continue
        
        full_spec_1 = data[:, i, j]
        
        if np.isnan(full_spec_1).any():
            print("NAN IN SPECTRUM")
            if multicomponent:
                fluxes_t.append(np.nan)
                fluxes_amp_1.append(np.nan)
                line_center_1.append(np.nan)
                fwhm_1.append(np.nan)
                fluxes_amp_2.append(np.nan)
                line_center_2.append(np.nan)
                fwhm_2.append(np.nan)
                if type == "triple":
                    fluxes_amp_3.append(np.nan)
                    line_center_3.append(np.nan)
                    fwhm_3.append(np.nan)
                continue
        
        if multicomponent:
            #print(full_spec_1)
            loc = LocalFit(wavelengths, full_spec_1, line_dict[name][2], line_dict[name][3], name)
            npoly = 1
            if type == "triple":
                loc.multicomponent_triple_fit(npoly=npoly)
            if type == "double":
                loc.multicomponent_double_fit(npoly=npoly)
        
            #loc.neon3_fit()
            #if i == 15 and j == 30:
            #    loc.render_fit(path="./test.pdf")
            #print(loc.popt)

            print(f"SPAXEL: {i}, {j}")
            
            #print(loc.line_strength.value)
            if loc.line_strength == 0:
                fluxes_t.append(np.nan)
                fluxes_amp_1.append(np.nan)
                line_center_1.append(np.nan)
                fwhm_1.append(np.nan)
                fluxes_amp_2.append(np.nan)
                line_center_2.append(np.nan)
                fwhm_2.append(np.nan)
                if type == "triple":
                    fluxes_amp_3.append(np.nan)
                    line_center_3.append(np.nan)
                    fwhm_3.append(np.nan)
                fluxes[i, j] = 0
            else:
                fluxes_t.append(loc.line_strength.value)
                fluxes_amp_1.append(loc.popt[npoly + 1])
                line_center_1.append(loc.popt[npoly + 2])
                fwhm_1.append(loc.popt[npoly + 3] * 2.355)
                fluxes_amp_2.append(loc.popt[npoly + 3 + 1])
                line_center_2.append(loc.popt[npoly + 3 + 2])
                fwhm_2.append(loc.popt[npoly + 3 + 3] * 2.355)
                if type == "triple":
                    fluxes_amp_3.append(loc.popt[npoly + 6 + 1])
                    line_center_3.append(loc.popt[npoly + 6 + 2])
                    fwhm_3.append(loc.popt[npoly + 6 + 3] * 2.355)
                fluxes[i, j] = loc.line_strength.value
        
        else:
            loc = LocalFit(wavelengths, full_spec_1, line_dict[name][2], line_dict[name][3], name)
            loc.main_fit(npoly=1, ngauss=1, spaxel_fit=True)
            print(f"SPAXEL: {i}, {j}")
            
            if loc.line_strength == 0:
                fluxes_t.append(np.nan)
                line_center.append(np.nan)
                fwhm.append(np.nan)
                fluxes[i, j] = 0
            else:
                fluxes_t.append(loc.line_strength.value)
                line_center.append(loc.popt[2] - line_dict[name][3])
                fwhm.append(loc.popt[3] * 2.355)
                fluxes[i, j] = loc.line_strength.value

if multicomponent:
    if type == "double":
        fitparams = Table([x_pix, y_pix, fluxes_t, fluxes_amp_1, line_center_1, fwhm_1, fluxes_amp_2, line_center_2, fwhm_2], names=("XPIX", "YPIX", "FLUXT", "AMP1", "DLINE1", "FWHM1", "AMP2", "DLINE2", "FWHM2"))
        fitparams.write(f"./../diagnostic_plots/spaxel_maps/{core}_{name}_{loc.npoly}_{loc.ngauss}_multi_double.dat", format="ipac", overwrite=True)
    if type == "triple":
        fitparams = Table([x_pix, y_pix, fluxes_t, fluxes_amp_1, line_center_1, fwhm_1, fluxes_amp_2, line_center_2, fwhm_2, fluxes_amp_3, line_center_3, fwhm_3], names=("XPIX", "YPIX", "FLUXT", "AMP1", "DLINE1", "FWHM1", "AMP2", "DLINE2", "FWHM2", "AMP3", "DLINE3", "FWHM3"))
        fitparams.write(f"./../diagnostic_plots/spaxel_maps/{core}_{name}_{loc.npoly}_{loc.ngauss}_multi_triple.dat", format="ipac", overwrite=True)
        
else:
    fitparams = Table([x_pix, y_pix, fluxes_t, line_center, fwhm], names=("XPIX", "YPIX", "FLUX", "LINEC", "FWHM"))
    fitparams.write(f"./../diagnostic_plots/spaxel_maps/{core}_{name}_{loc.npoly}_{loc.ngauss}.dat", format="ipac", overwrite=True)

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