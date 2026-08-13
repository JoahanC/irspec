import glob 
from irspec.datacube import Datacube
from irspec.cubespec import CubeSpec
import astropy.units as u
import matplotlib.pyplot as plt
from irspec.emission_io import read_line_params3

test_path = glob.glob(f"./irspec/north_data/ch1-short_s3d.fits")[0]
testcube = Datacube(test_path, redshift=0.044601, verbose=True)
spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                    "ir23128s_single2S_param.txt", "./../../creta_extractions/single2S/",
                    "./../../cafe_output/single2S_SB_and_AGN/", 0.044601, testcube, mode="SB_and_AGN")
lines = ["[FeII]", "[ArII]", "[SIII]", "[ArIII]", "[NeII]", "[SIV]", "[NeIII]"]
lines_ip = [16.19, 27.63, 34.79, 40.74, 40.96, 47.22, 63.45]
coronal_lines = ["[OIV]", "[NeV]", "[NeV]_14"]
coronal_lines_ip = [77.41, 126.21, 126.21]
all_lines = ["[FeII]", "[ArII]", "[SIII]", "[ArIII]", "[NeII]", "[SIV]", "[NeIII]", "[OIV]", "[NeV]", "[NeV]_14"]
all_lines = ["[ArII]", "[SIII]", "[ArIII]", "[NeII]", "[SIV]", "[NeIII]", "[OIV]", "[NeV]", "[NeV]_14"]

all_lines_ip = [16.19, 27.63, 34.79, 40.74, 40.96, 47.22, 63.45, 77.41, 126.21, 126.21]
all_lines_ip = [27.63, 34.79, 40.74, 40.96, 47.22, 63.45, 77.41, 126.21, 126.21]

molhy_lines = ["[H_2_S_1]", "[H_2_S_2]", "[H_2_S_3]", "[H_2_S_5]"]
line_dict = read_line_params3()

colors = ["red", "blue", "green", "orange", "magenta", "brown", "purple", "cyan", "grey", "orangered", "indigo"]

normal = True
coronal = False 
combined = True
mol2 = False
ratios = True
if normal:
    fig, ax = plt.subplots()
    for idx, line in enumerate(lines):
        values = spec_obj.refit_lines([lines[idx]])
        amps = values[0]
        rel_vels = values[1]
        vel_disp = values[2]
        narrow_idx = 0
        broad_idx = 1
        if len(rel_vels) == 2:
            if vel_disp[0].value < vel_disp[1].value:
                narrow_idx = 0
                broad_idx = 1
            else:
                broad_idx = 0
                narrow_idx = 1
            ax.scatter(lines_ip[idx], vel_disp[narrow_idx].value, label=line, marker="^", c=colors[idx], s=100)
            ax.scatter(lines_ip[idx], vel_disp[broad_idx].value, label=line, marker="o", c=colors[idx], s=100)
        else:
            ax.scatter(lines_ip[idx], vel_disp[0].value, label=line, marker="^", c=colors[idx], s=100)
    ax.legend()
    ax.set_xlabel("Ionization Potential [eV]")
    ax.set_ylabel(r"$\sigma$ [km/s]")
    plt.show()

if coronal:
    fig, ax = plt.subplots()
    for idx, line in enumerate(coronal_lines):
        values = spec_obj.refit_lines([coronal_lines[idx]])
        amps = values[0]
        rel_vels = values[1]
        vel_disp = values[2]
        narrow_idx = 0
        broad_idx = 1
        if len(rel_vels) == 2:
            if vel_disp[0].value < vel_disp[1].value:
                narrow_idx = 0
                broad_idx = 1
            else:
                broad_idx = 0
                narrow_idx = 1
            ax.scatter(coronal_lines_ip[idx], vel_disp[narrow_idx].value, label=line, marker="^", c=colors[idx], s=100)
            ax.scatter(coronal_lines_ip[idx], vel_disp[broad_idx].value, label=line, marker="o", c=colors[idx], s=100)
        else:
            ax.scatter(coronal_lines_ip[idx], vel_disp[0].value, label=line, marker="^", c=colors[idx], s=100)
    ax.legend()
    ax.set_xlabel("Ionization Potential [eV]")
    ax.set_ylabel(r"$\sigma$ [km/s]")
    plt.show()

if combined:
    fig, ax = plt.subplots()
    for idx, line in enumerate(all_lines):
        values = spec_obj.refit_lines([all_lines[idx]])
        amps = values[0]
        rel_vels = values[1]
        vel_disp = values[2]
        narrow_idx = 0
        broad_idx = 1
        if len(rel_vels) == 2:
            if vel_disp[0].value < vel_disp[1].value:
                narrow_idx = 0
                broad_idx = 1
            else:
                broad_idx = 0
                narrow_idx = 1
            ax.scatter(all_lines_ip[idx], vel_disp[narrow_idx].value, label=line, marker="^", c=colors[idx], s=100)
            ax.scatter(all_lines_ip[idx], vel_disp[broad_idx].value, label=line, marker="o", c=colors[idx], s=100)
        else:
            if len(vel_disp) != 0:
                ax.scatter(all_lines_ip[idx], vel_disp[0].value, label=line, marker="^", c=colors[idx], s=100)
    ax.legend()
    ax.set_xlabel("Ionization Potential [eV]")
    ax.set_ylabel(r"$\sigma$ [km/s]")
    plt.show()

if mol2:
    for idx, line in enumerate(molhy_lines):
        values = spec_obj.refit_lines([molhy_lines[idx]])
        amps = values[0]
        rel_vels = values[1]
        vel_disp = values[2]
        narrow_idx = 0
        broad_idx = 1

ratio_top = ["[NeV]_14", "[SIV]", "[NeIII]", "[OIV]", "[NeIII]"]
ratio_bot = ["[NeII]", "[NeII]", "[NeII]", "[NeII]", "[NeV]_14"]
if ratios:
    for idx, line in enumerate(ratio_top):
        fluxes_top = spec_obj.refit_lines([ratio_top[idx]])[0]
        fluxes_bot = spec_obj.refit_lines([ratio_bot[idx]])[0]
        tot_flux_top = 0
        tot_flux_bot = 0
        for flux in fluxes_top:
            tot_flux_top += flux.value
        for flux in fluxes_bot:
            tot_flux_bot += flux.value
        print(f"{ratio_top[idx]}/{ratio_bot[idx]}")
        print(tot_flux_top / tot_flux_bot)