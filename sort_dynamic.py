""" 
This script sorts the fitted values returned from the dynamic multicomponent fitting.
"""

from astropy.io import ascii
from astropy.table import Table 
import numpy as np


def sortdata(sorting_method, name):
    data = ascii.read(f"./../diagnostic_plots/dynamic_multicomponent/{name}/SNRCUT_fit.dat", format="ipac")  

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
        g1_amp = data["G1AMP"][idx]
        g2_amp = data["G2AMP"][idx]
        g3_amp = data["G3AMP"][idx]
        g1_vel = data["G1CEN"][idx]
        g2_vel = data["G2CEN"][idx]
        g3_vel = data["G3CEN"][idx]
        g1_sigma = data["G1SIGMA"][idx]
        g2_sigma = data["G2SIGMA"][idx]
        g3_sigma = data["G3SIGMA"][idx]
        
        if np.isnan(g3_amp):
            if np.isnan(g2_amp):
                if np.isnan(g1_amp):
                    g1_amps.append(np.nan)
                    g2_amps.append(np.nan)
                    g3_amps.append(np.nan)
                    g1_vels.append(np.nan)
                    g2_vels.append(np.nan)
                    g3_vels.append(np.nan)
                    g1_sigmas.append(np.nan)
                    g2_sigmas.append(np.nan)
                    g3_sigmas.append(np.nan)
                else:
                    g1_amps.append(g1_amp)
                    g2_amps.append(np.nan)
                    g3_amps.append(np.nan)
                    g1_vels.append(g1_vel)
                    g2_vels.append(np.nan)
                    g3_vels.append(np.nan)
                    g1_sigmas.append(g1_sigma)
                    g2_sigmas.append(np.nan)
                    g3_sigmas.append(np.nan)
            else:
                if sorting_method == "SIGMA":
                    if g1_sigma > g2_sigma:
                        g1_amps.append(g1_amp)
                        g2_amps.append(g2_amp)
                        g3_amps.append(np.nan)
                        g1_vels.append(g1_vel)
                        g2_vels.append(g2_vel)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(g1_sigma)
                        g2_sigmas.append(g2_sigma)
                        g3_sigmas.append(np.nan)
                    else:
                        g1_amps.append(g2_amp)
                        g2_amps.append(g1_amp)
                        g3_amps.append(np.nan)
                        g1_vels.append(g2_vel)
                        g2_vels.append(g1_vel)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(g2_sigma)
                        g2_sigmas.append(g1_sigma)
                        g3_sigmas.append(np.nan)
                
                if sorting_method == "CEN":
                    if g1_vel < g2_vel:
                        g1_amps.append(g1_amp)
                        g2_amps.append(g2_amp)
                        g3_amps.append(np.nan)
                        g1_vels.append(g1_vel)
                        g2_vels.append(g2_vel)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(g1_sigma)
                        g2_sigmas.append(g2_sigma)
                        g3_sigmas.append(np.nan)
                    else:
                        g1_amps.append(g2_amp)
                        g2_amps.append(g1_amp)
                        g3_amps.append(np.nan)
                        g1_vels.append(g2_vel)
                        g2_vels.append(g1_vel)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(g2_sigma)
                        g2_sigmas.append(g1_sigma)
                        g3_sigmas.append(np.nan)
                
                if sorting_method == "AMP":
                    if g1_amp > g2_amp:
                        g1_amps.append(g1_amp)
                        g2_amps.append(g2_amp)
                        g3_amps.append(np.nan)
                        g1_vels.append(g1_vel)
                        g2_vels.append(g2_vel)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(g1_sigma)
                        g2_sigmas.append(g2_sigma)
                        g3_sigmas.append(np.nan)
                    else:
                        g1_amps.append(g2_amp)
                        g2_amps.append(g1_amp)
                        g3_amps.append(np.nan)
                        g1_vels.append(g2_vel)
                        g2_vels.append(g1_vel)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(g2_sigma)
                        g2_sigmas.append(g1_sigma)
                        g3_sigmas.append(np.nan)
        
        # ALL IN PLAY
        else:
            if sorting_method == "SIGMA":
                dict_keys = [g1_sigma, g2_sigma, g3_sigma]
                sort_dict = {g1_sigma: [g1_amp, g1_vel], g2_sigma: [g2_amp, g2_vel], g3_sigma: [g3_amp, g3_vel]}
                dict_keys.sort(reverse=True)
                g1_amps.append(sort_dict[dict_keys[0]][0])
                g2_amps.append(sort_dict[dict_keys[1]][0])
                g3_amps.append(sort_dict[dict_keys[2]][0])
                g1_vels.append(sort_dict[dict_keys[0]][1])
                g2_vels.append(sort_dict[dict_keys[1]][1])
                g3_vels.append(sort_dict[dict_keys[2]][1])
                g1_sigmas.append(dict_keys[0])
                g2_sigmas.append(dict_keys[1])
                g3_sigmas.append(dict_keys[2])
            
            if sorting_method == "CEN":
                dict_keys = [g1_vel, g2_vel, g3_vel]
                sort_dict = {g1_vel: [g1_amp, g1_sigma], g2_vel: [g2_amp, g2_sigma], g3_vel: [g3_amp, g3_sigma]}
                dict_keys.sort(reverse=True)
                g1_amps.append(sort_dict[dict_keys[0]][0])
                g2_amps.append(sort_dict[dict_keys[1]][0])
                g3_amps.append(sort_dict[dict_keys[2]][0])
                g1_vels.append(g1_vel)
                g2_vels.append(g2_vel)
                g3_vels.append(g3_vel)
                g1_sigmas.append(sort_dict[dict_keys[0]][1])
                g2_sigmas.append(sort_dict[dict_keys[1]][1])
                g3_sigmas.append(sort_dict[dict_keys[2]][1])
            
            if sorting_method == "AMP":
                dict_keys = [g1_amp, g2_amp, g3_amp]
                sort_dict = {g1_amp: [g1_vel, g1_sigma], g2_amp: [g2_vel, g2_sigma], g3_amp: [g3_vel, g3_sigma]}
                dict_keys.sort(reverse=True)
                g1_amps.append(g1_amp)
                g2_amps.append(g2_amp)
                g3_amps.append(g3_amp)
                g1_vels.append(sort_dict[dict_keys[0]][0])
                g2_vels.append(sort_dict[dict_keys[1]][0])
                g3_vels.append(sort_dict[dict_keys[2]][0])
                g1_sigmas.append(sort_dict[dict_keys[0]][1])
                g2_sigmas.append(sort_dict[dict_keys[1]][1])
                g3_sigmas.append(sort_dict[dict_keys[2]][1])
        
    best_fit_values = [data["XPIX"], data["YPIX"], g1_amps, g1_vels, g1_sigmas, g2_amps, g2_vels, g2_sigmas, g3_amps, g3_vels, g3_sigmas]
    fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA"))
    fitparams.write(f"./../diagnostic_plots/dynamic_multicomponent/{name}/{sorting_method}_fit.dat", format="ipac", overwrite=True)

sorting_method = "CEN"
name = "[NeV]_14"
sortdata(sorting_method, name)