
from astropy.io import ascii
from astropy.table import Table 
import numpy as np

line_name = "[NeV]_14"
line_center = 14.3217
data = ascii.read(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/AMP_fit.dat", format="ipac")  
print(data)

g1_amps = []
g2_amps = []
g3_amps = []
g1_vels = []
g2_vels = []
g3_vels = []
g1_sigmas = []
g2_sigmas = []
g3_sigmas = []

for idx, _ in enumerate(data["XPIX"]):
    # both nan
    if np.isnan(data["G1CEN"][idx]) and np.isnan(data["G2CEN"][idx]):
        g1_amps.append(data["G1AMP"][idx])
        g2_amps.append(data["G2AMP"][idx])
        g1_vels.append(data["G1CEN"][idx])
        g2_vels.append(data["G2CEN"][idx])
        g1_sigmas.append(data["G1SIGMA"][idx])
        g2_sigmas.append(data["G2SIGMA"][idx])

    if np.abs(data["G1CEN"][idx] - line_center) > np.abs(data["G2CEN"][idx] - line_center):
        g1_amps.append(data["G1AMP"][idx])
        g2_amps.append(data["G2AMP"][idx])
        g1_vels.append(data["G1CEN"][idx])
        g2_vels.append(data["G2CEN"][idx])
        g1_sigmas.append(data["G1SIGMA"][idx])
        g2_sigmas.append(data["G2SIGMA"][idx])
    if np.abs(data["G1CEN"][idx] - line_center) <= np.abs(data["G2CEN"][idx] - line_center):
        g1_amps.append(data["G2AMP"][idx])
        g2_amps.append(data["G1AMP"][idx])
        g1_vels.append(data["G2CEN"][idx])
        g2_vels.append(data["G1CEN"][idx])
        g1_sigmas.append(data["G2SIGMA"][idx])
        g2_sigmas.append(data["G1SIGMA"][idx])
    # Check if only 2 element is nan
    if not np.isnan(data["G1CEN"][idx]) and np.isnan(data["G2CEN"][idx]):
        g1_amps.append(data["G1AMP"][idx])
        g2_amps.append(data["G2AMP"][idx])
        g1_vels.append(data["G1CEN"][idx])
        g2_vels.append(data["G2CEN"][idx])
        g1_sigmas.append(data["G1SIGMA"][idx])
        g2_sigmas.append(data["G2SIGMA"][idx])
        
    g3_amps.append(data["G3AMP"][idx])
    g3_vels.append(data["G3CEN"][idx])
    g3_sigmas.append(data["G3SIGMA"][idx])

print(len(data["XPIX"]), len(g1_amps), len(g2_amps), len(g3_amps), len(g1_vels), len(g2_vels), len(g3_vels) , len(g1_sigmas), len(g2_sigmas), len(g3_sigmas))
best_fit_values = [data["XPIX"], data["YPIX"], g1_amps, g1_vels, g1_sigmas, g2_amps, g2_vels, g2_sigmas, g3_amps, g3_vels, g3_sigmas]
fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA"))
fitparams.write(f"./../diagnostic_plots/dynamic_multicomponent/{line_name}/AMP_CEN_fit.dat", format="ipac", overwrite=True)