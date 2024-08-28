import matplotlib.pyplot as plt 
import numpy as np
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
import os
import glob
from astroquery.ipac.irsa import Irsa
from astroquery import sha
from astropy.table import Table



"""filepath = "./../archival_data/spitzer_spectra/ir23128"
spitzer_files = os.listdir(filepath)
for filename in spitzer_files:
    path = os.path.join(filepath, filename)
    fits.open(path, ignore_missing_simple=True)
    #
    #data = fits.getdata(path)
    #hdu = fits.open(path, ignore_missing_simple=True, ext=1)[1]
    #print(data)"""

"""catlist = Irsa.list_catalogs()
for cat in catlist:
    #if "spit" in cat:
    #    print(cat)
    if cat == "spitzer.goals_irs_spec":
        print(cat, catlist[cat])
#"""

"""coord = coord.SkyCoord(ra=348.94466, dec=-59.0530, unit=(u.degree, u.degree))
table = Irsa.query_region(coordinates=coord, catalog="spitzer.goals_irs_spec", spatial="Cone",
                          radius=2 * u.arcmin)
print(table.keys())
print(table["spectrum_filename"])
print(table["spectrum_download_u"])

sha.get_file("https://" + table["spectrum_download_u"][0])"""

spitzer_files = glob.glob("./../archival_data/spitzer_spectra/ir23128/*.tbl")
#print(spitzer_files)
testfile = spitzer_files[0]



"""fig, ax = plt.subplots()
for idx, file in enumerate(spitzer_files):
    t = Table.read(file, format='ipac')
    ax.plot(t["wavelength"], t["flux_density"], label=idx)
    ax.set_yscale("log")
ax.legend()
plt.show()"""

def remove_outliers(fluxes):
    
    corrected_fluxes = np.copy(fluxes)
    
    for idx in range(len(fluxes) - 1):
        if idx == 0:
            continue
        past_flux = fluxes[idx - 1]
        curr_flux = fluxes[idx]
        fut_flux = fluxes[idx + 1]
        if curr_flux < (past_flux + fut_flux) / 2 * 0.5:
            corrected_fluxes[idx] = (past_flux + fut_flux) / 2
    
    return corrected_fluxes

fig, ax = plt.subplots()
t = Table.read(testfile, format='ipac')
ax.plot(t["wavelength"], t["flux_density"])
corrected_flux = remove_outliers(t["flux_density"])
ax.plot(t["wavelength"], corrected_flux)
for idx, _ in enumerate(t["wavelength"]):
    ax.scatter(t["wavelength"][idx], t["flux_density"][idx])
    #ax.text(t["wavelength"][idx], t["flux_density"][idx] * 1.1, idx, bbox=dict(boxstyle="round",
    #                ec=(1., 1, 1),
    #                fc=(0., 0, 1),
    #                ))
ax.set_yscale("log")
ax.legend()
plt.show()




def calculateStitchRatios(realData, aperture_correction, idx, grid):

    allRatio = []
    for i in range(len(realData)-1):
        
        chA = realData[i]
        chB = realData[i+1]
        
        if aperture_correction: 
            #print('aderfia edw eimaste complettt')
            if grid:
                chA_data = np.array(chA.spectrum_PSF_corrected) # not tested
                chB_data = np.array(chB.spectrum_PSF_corrected) # not tested
            else:
                chA_data = np.array(chA.spectrum_PSF_corrected)[idx,:]
                chB_data = np.array(chB.spectrum_PSF_corrected)[idx,:]
                
        else:   
            if grid:
                chA_data = np.array(chA.corrected_spectrum)[idx,:]
                chB_data = np.array(chB.corrected_spectrum)[idx,:]
            else:    
                chA_data = np.array(chA.corrected_spectrum)[:,idx]
                chB_data = np.array(chB.corrected_spectrum)[:,idx]
                
        chA_ls = chA.ls
        chB_ls = chB.ls
        delta1 = chA.CDELT3
        delta2 = chB.CDELT3
        chA_start = np.where(np.array(chA_ls) >= chB_ls[0] - delta2/2)[0]
        chB_stop =  np.where(np.array(chB_ls) <= chA_ls[len(chA_ls)-1] + delta1/2)[0]
        
        if len(chA_start) == 0 or len(chB_stop) == 0:
            ratio = 1.
        
        else:
            chA_overlapping_data = chA_data[chA_start[0]:]
            chA_over_ls= chA_ls[chA_start[0]:]
            chB_overlapping_data = chB_data[:chB_stop[-1]+1]
            
            chA_median = np.nanmedian(chA_overlapping_data)
            chB_median = np.nanmedian(chB_overlapping_data)
            
            if np.isnan(chA_median) or np.isnan(chB_median):
                ratio = 1.
            else:
                ratio = chB_median / chA_median 
                #chA_mean = np.mean(chA_overlapping_data)
                #chB_mean = np.mean(chB_overlapping_data)
                #ratio = chB_mean / chA_mean 
        
        allRatio.append(ratio)

    allRatio.append(1.)

    return allRatio