""" 
This file writes a basic data file to aggregate published GOALS sample
data
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from plotparams import PlotParams


# Load in archival data files
st11_path = "./../archival_data/J_ApJS_206_1/table1.dat"
st21_path = "./../archival_data/J_ApJ_777_156/table1.dat"
st22_path = "./../archival_data/J_ApJ_777_156/table2.dat"
st23_path = "./../archival_data/J_ApJ_777_156/table3.dat"
st31_path = "./../archival_data/J_ApJ_790_124/table1.dat"
st32_path = "./../archival_data/J_ApJ_790_124/table1.dat"




goals_sample = {}

names = []
merger_stage = []
with open(st11_path, 'r') as file:
    data11 = file.readlines()
    for line in data11:
        names.append(line[:17].strip())
        merger_stage.append(line[119])
for idx, name in enumerate(names):
    if name not in goals_sample:
        goals_sample[name] = []

st1names = []
with open(st21_path, 'r') as file:
    data21 = file.readlines()
    for line in data21:
        name = line[:17].strip()
        st1names.append(name)
        if name not in goals_sample:
            goals_sample[name] = []
        
        s4flux = line[58:63].strip()
        if s4flux != "":
            s4flux = float(s4flux)
        if s4flux == "":
            s4flux = np.nan
        if s4flux == np.nan:
            s4_flux_err = np.nan 
        if s4flux != np.nan:
            s4_flux_err = line[65:69].strip()
            s4_flux_err = float(s4_flux_err) 
        
        ne2flux = line[72:77].strip()
        if ne2flux != "":
            ne2flux = float(ne2flux)
        if ne2flux == "":
            ne2flux = np.nan 
        """if ne2flux == np.nan:
            ne2_flux_err = np.nan 
        if ne2flux != np.nan:
            ne2_flux_err = line[79:83].strip()
            ne2_flux_err = float(ne2_flux_err)"""
        
        ne5flux = line[86:91].strip()
        if ne5flux != "":
            ne5flux = float(ne5flux)
        if ne5flux == "":
            ne5flux = np.nan 
        """if ne5flux == np.nan:
            ne5_flux_err = np.nan 
        if ne5flux != np.nan:
            ne5_flux_err = line[93:97].strip()
            ne5_flux_err = float(ne5_flux_err)"""
        
        ne3flux = line[100:106].strip()
        if ne3flux != "":
            ne3flux = float(ne3flux)
        if ne3flux == "":
            ne3flux = np.nan
        """if ne3flux == np.nan:
            ne3_flux_err = np.nan 
        if ne3flux != np.nan:
            ne3_flux_err = line[108:112].strip()
            ne3_flux_err = float(ne3_flux_err) """
        
        goals_sample[name].append(s4flux)
        #goals_sample[name].append(s4_flux_err)
        goals_sample[name].append(ne2flux)
        #goals_sample[name].append(ne2_flux_err)
        goals_sample[name].append(ne5flux)
        #goals_sample[name].append(ne5_flux_err)
        goals_sample[name].append(ne3flux)
        #goals_sample[name].append(ne3_flux_err)

for idx, name in enumerate(goals_sample):
    if name=="IRASF15250+3608":
        print(idx)
    if name=="ESO286-IG019":
        print(idx)

ne2_fluxes = []
ne3_fluxes = []
ne5_fluxes = []
s4_fluxes = []
for goal in goals_sample:
    if len(goals_sample[goal]) != 4:
        for _ in range(4):
            goals_sample[goal].append(np.nan)
    s4_fluxes.append(goals_sample[goal][0])
    ne2_fluxes.append(goals_sample[goal][2])
    ne5_fluxes.append(goals_sample[goal][3])
    ne3_fluxes.append(goals_sample[goal][4])

ne2_fluxes = np.array(ne2_fluxes)
ne3_fluxes = np.array(ne3_fluxes)
ne5_fluxes = np.array(ne5_fluxes)
s4_fluxes = np.array(s4_fluxes)

ne3_ne2_log_ratio = np.log10(ne3_fluxes / ne2_fluxes)
s4_ne2_log_ratio = np.log10(s4_fluxes / ne2_fluxes)

plt.style.use('dark_background')
fig, ax = plt.subplots()

fig.set_size_inches(8, 8)
ax.scatter(ne3_ne2_log_ratio, s4_ne2_log_ratio, s=10, color="white")
ax.set_xlabel("log([NeIII]/Ne[II])")
ax.set_ylabel("log([SIV]/[NeII])")
plt.savefig("./../diagnostic_plots/s4_ne3_ne2_ratios.pdf", dpi=1000)
plt.close()

ir15250_line_path = "./../cafe_output/IR15250/IR15250_SingleExt_r03as/IR15250_SingleExt_r03as_linetable.ecsv"
ir23128n_line_path = "./../cafe_output/IR23128-N/IR23128-N_SingleExt_r03as/IR23128-N_SingleExt_r03as_linetable.ecsv"
ir23128s_line_path = "./../cafe_output/IR23128-S/IR23128-S_SingleExt_r07as/IR23128-S_SingleExt_r07as_linetable.ecsv"

ir15250_line_table = Table.read(ir15250_line_path, format='ascii.ecsv')
ir23128n_line_table = Table.read(ir23128n_line_path, format='ascii.ecsv')
ir23128s_line_table = Table.read(ir23128s_line_path, format='ascii.ecsv')

# IR15250
ir15250_narrow_s4_idx = np.argwhere(ir15250_line_table["line_name"] == "SIV_105105N")[0][0]
ir15250_broad_s4_idx = np.argwhere(ir15250_line_table["line_name"] == "SIV_105105B")[0][0]
ir15250_ne2_idx = np.argwhere(ir15250_line_table["line_name"] == "NeII_128136N")[0][0]
ir15250_narrow_ne3_idx = np.argwhere(ir15250_line_table["line_name"] == "NeIII_155551N")[0][0]
ir15250_broad_ne3_idx = np.argwhere(ir15250_line_table["line_name"] == "NeIII_155551B")[0][0]

ir15250_s4_flux = ir15250_line_table["line_strength_int"][ir15250_broad_s4_idx] + ir15250_line_table["line_strength_int"][ir15250_narrow_s4_idx]
ir15250_s4_flux_err = np.sqrt(ir15250_line_table["line_strength_int_unc"][ir15250_broad_s4_idx] ** 2 + ir15250_line_table["line_strength_int_unc"][ir15250_narrow_s4_idx] ** 2)
ir15250_ne2_flux = ir15250_line_table["line_strength_int"][ir15250_ne2_idx]
ir15250_ne2_flux_err = ir15250_line_table["line_strength_int_unc"][ir15250_ne2_idx]
ir15250_ne3_flux = ir15250_line_table["line_strength_int"][ir15250_broad_ne3_idx] + ir15250_line_table["line_strength_int"][ir15250_narrow_ne3_idx]
ir15250_ne3_flux_err = np.sqrt(ir15250_line_table["line_strength_int_unc"][ir15250_broad_ne3_idx] ** 2 + ir15250_line_table["line_strength_int_unc"][ir15250_narrow_ne3_idx] ** 2)

# IR23128N
ir23128n_narrow_s4_idx = np.argwhere(ir23128n_line_table["line_name"] == "SIV_105105N")[0][0]
ir23128n_broad_s4_idx = np.argwhere(ir23128n_line_table["line_name"] == "SIV_105105B")[0][0]
ir23128n_ne2_idx = np.argwhere(ir23128n_line_table["line_name"] == "NeII_128136N")[0][0]
ir23128n_narrow_ne3_idx = np.argwhere(ir23128n_line_table["line_name"] == "NeIII_155551N")[0][0]
ir23128n_broad_ne3_idx = np.argwhere(ir23128n_line_table["line_name"] == "NeIII_155551B")[0][0]

ir23128n_s4_flux = ir23128n_line_table["line_strength_int"][ir23128n_broad_s4_idx] + ir23128n_line_table["line_strength_int"][ir23128n_narrow_s4_idx]
ir23128n_s4_flux_err = np.sqrt(ir23128n_line_table["line_strength_int_unc"][ir23128n_broad_s4_idx] ** 2 + ir23128n_line_table["line_strength_int_unc"][ir23128n_narrow_s4_idx] ** 2)
ir23128n_ne2_flux = ir23128n_line_table["line_strength_int"][ir23128n_ne2_idx]
ir23128n_ne2_flux_err = ir23128n_line_table["line_strength_int_unc"][ir23128n_ne2_idx]
ir23128n_ne3_flux = ir23128n_line_table["line_strength_int"][ir23128n_broad_ne3_idx] + ir23128n_line_table["line_strength_int"][ir23128n_narrow_ne3_idx]
ir23128n_ne3_flux_err = np.sqrt(ir23128n_line_table["line_strength_int_unc"][ir23128n_broad_ne3_idx] ** 2 + ir23128n_line_table["line_strength_int_unc"][ir23128n_narrow_ne3_idx] ** 2)

# IR23128S
ir23128s_narrow_s4_idx = np.argwhere(ir23128s_line_table["line_name"] == "SIV_105105N")[0][0]
ir23128s_broad_s4_idx = np.argwhere(ir23128s_line_table["line_name"] == "SIV_105105B")[0][0]
ir23128s_ne2_idx = np.argwhere(ir23128s_line_table["line_name"] == "NeII_128136N")[0][0]
ir23128s_narrow_ne3_idx = np.argwhere(ir23128s_line_table["line_name"] == "NeIII_155551N")[0][0]
ir23128s_broad_ne3_idx = np.argwhere(ir23128s_line_table["line_name"] == "NeIII_155551B")[0][0]

ir23128s_s4_flux = ir23128s_line_table["line_strength_int"][ir23128s_broad_s4_idx] + ir23128s_line_table["line_strength_int"][ir23128s_narrow_s4_idx]
ir23128s_s4_flux_err = np.sqrt(ir23128s_line_table["line_strength_int_unc"][ir23128s_broad_s4_idx] ** 2 + ir23128s_line_table["line_strength_int_unc"][ir23128s_narrow_s4_idx] ** 2)
ir23128s_ne2_flux = ir23128s_line_table["line_strength_int"][ir23128s_ne2_idx]
ir23128s_ne2_flux_err = ir23128s_line_table["line_strength_int_unc"][ir23128s_ne2_idx]
ir23128s_ne3_flux = ir23128s_line_table["line_strength_int"][ir23128s_broad_ne3_idx] + ir23128s_line_table["line_strength_int"][ir23128s_narrow_ne3_idx]
ir23128s_ne3_flux_err = np.sqrt(ir23128s_line_table["line_strength_int_unc"][ir23128s_broad_ne3_idx] ** 2 + ir23128s_line_table["line_strength_int_unc"][ir23128s_narrow_ne3_idx] ** 2)

new_s4_fluxes = np.array([ir15250_s4_flux, ir23128n_s4_flux, ir23128s_s4_flux])
new_s4_flux_errs = np.array([ir15250_s4_flux_err, ir23128n_s4_flux_err, ir23128s_s4_flux_err])
new_ne2_fluxes = np.array([ir15250_ne2_flux, ir23128n_ne2_flux, ir23128s_ne2_flux])
new_ne2_flux_errs = np.array([ir15250_ne2_flux_err, ir23128n_ne2_flux_err, ir23128s_ne2_flux_err])
new_ne3_fluxes = np.array([ir15250_ne3_flux, ir23128n_ne3_flux, ir23128s_ne3_flux])
new_ne3_flux_errs = np.array([ir15250_ne3_flux_err, ir23128n_ne3_flux_err, ir23128s_ne3_flux_err])

new_ne3_ne2_log_ratio = np.log10(new_ne3_fluxes / new_ne2_fluxes)
new_ne3_ne2_log_ratio_err = np.log10(np.sqrt(new_ne3_flux_errs ** 2 + new_ne2_flux_errs ** 2))
new_s4_ne2_log_ratio = np.log10(new_s4_fluxes / new_ne2_fluxes)
new_s4_ne2_log_ratio_err = np.log10(np.sqrt(new_s4_flux_errs ** 2 + new_ne2_flux_errs ** 2))

plt.style.use('dark_background')
fig, ax = plt.subplots()

print(ne3_fluxes[252])

fig.set_size_inches(9, 8)
ax.scatter(ne3_ne2_log_ratio, s4_ne2_log_ratio, s=50, alpha=0.5, color="white", label="Inami et al. (2013)")
ax.scatter(ne3_ne2_log_ratio[162], s4_ne2_log_ratio[162], s=300, marker="D", color="lime", label="Inami - IR15250")
ax.scatter(ne3_ne2_log_ratio[252], s4_ne2_log_ratio[252], s=300, marker="D", color="gold", label="Inami - IR23128")
ax.scatter(new_ne3_ne2_log_ratio[0], new_s4_ne2_log_ratio[0], s=300, marker="*", color="lime", label="JWST - IR15250")
ax.scatter(new_ne3_ne2_log_ratio[1], new_s4_ne2_log_ratio[1], s=300, marker="*", color="fuchsia", label="JWST - IR23128-N")
ax.scatter(new_ne3_ne2_log_ratio[2], new_s4_ne2_log_ratio[2], s=300, marker="*", color="gold", label="JWST - IR23128-S")
ax.set_xlabel("log([NeIII]/Ne[II])", fontsize=24)
ax.set_ylabel("log([SIV]/[NeII])", fontsize=24)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)
plt.legend(prop={'size': 16})
plt.savefig("./../diagnostic_plots/new_s4_ne3_ne2_ratios.pdf", dpi=1000)
plt.close()

"""# neon line width
neon_dict = {"name": [], "ne2lw": [], "ne2lw_unc": [], "ne3lw": [], "ne3lw_unc": [], "ne5lw": [], "ne5lw_unc": []}
neon_line_path = "./../J_ApJ_777_156/table3.dat"
with open(neon_line_path, 'r') as file:
    for line in file.readlines():
        neon_dict["name"].append(line[:18].strip())
        ne2 = line[18:22].strip()
        if ne2 != '':
            neon_dict["ne2lw"].append(float(ne2))
        if ne2 == '':
            neon_dict["ne2lw"].append(np.nan)
        ne2_unc = line[22:25].strip()
        if ne2_unc != '':
            neon_dict["ne2lw_unc"].append(float(line[22:25].strip()))
        if ne2_unc == '':
            neon_dict["ne2lw_unc"].append(np.nan)
        ne3 = line[25:28].strip()
        if ne3 != '':
            neon_dict["ne3lw"].append(float(line[25:28].strip()))
        if ne3 == '':
            neon_dict["ne3lw"].append(np.nan)
        ne3_unc = line[29:32].strip()
        if ne3_unc != '':
            neon_dict["ne3lw_unc"].append(float(line[29:32].strip()))
        if ne3_unc == '':
            neon_dict["ne3lw_unc"].append(np.nan)
        ne5 = line[32:36].strip()
        if ne5 != '':
            neon_dict["ne5lw"].append(float(line[32:36].strip()))
        if ne5 == '':
            neon_dict["ne5lw"].append(np.nan)
        ne5_unc = line[37:].strip()
        if ne5_unc != '':
            neon_dict["ne5lw_unc"].append(float(line[37:].strip()))
        if ne5_unc == '':
            neon_dict["ne5lw_unc"].append(np.nan)"""


eqw_dict = {"name": [], "EQW6.2": [], "6.2f": [], "7.7f": [], "11.3f": []}
steirwalt_path = "./../archival_data/J_ApJ_790_124/table1.dat"
with open(steirwalt_path, 'r') as file:
    for line in file.readlines():
        eqw_dict["name"].append(line[:17].strip())
        eqw62 = line[24:28].strip()
        if eqw62 == "":
            eqw_dict["EQW6.2"].append(np.nan)
        if eqw62 != "":
            eqw_dict["EQW6.2"].append(float(eqw62))
        f62 = line[38:43].strip()
        if f62 == "":
            eqw_dict["6.2f"].append(np.nan)
        if f62 != "":
            eqw_dict["6.2f"].append(float(f62))
        f77 = line[51:57].strip()
        if f77 == "":
            eqw_dict["7.7f"].append(np.nan)
        if f77 != "":
            eqw_dict["7.7f"].append(float(f77))
        f113 = line[66:71].strip()
        if f113 == "":
            eqw_dict["11.3f"].append(np.nan)
        if f113 != "":
            eqw_dict["11.3f"].append(float(f113))


ir23128s_path = "./../cafe_output/IR23128-S/IR23128-S_SingleExt_r07as/IR23128-S_SingleExt_r07as_pahtable.ecsv"
ir23128n_path = "./../cafe_output/IR23128-N/IR23128-N_SingleExt_r03as/IR23128-N_SingleExt_r03as_pahtable.ecsv"
ir23128s_path = "./../cafe_output/IR15250/IR15250_SingleExt_r03as/IR15250_SingleExt_r03as_pahtable.ecsv"


fluxes_62_new = np.array([1.214630846166026e-15, 1.5151752618603605e-16, 3.3175474354322733e-16])
fluxes_77_new = np.array([3.893850351354496e-15, 4.94602697351627e-16, 6.921628356277496e-15])
fluxes_113_new = np.array([6.617338616963007e-16, 6.558926881852666e-17, 4.2045133453377543e-16])
eqw62_new = np.array([0.5981623983947937, 2.260544179417326, 0.31937122652915545])

ratio_113_77_new = fluxes_113_new / fluxes_77_new
ratio_62_77_new = fluxes_62_new / fluxes_77_new

eqw62 = np.array(eqw_dict["EQW6.2"])
fluxes_62 = np.array(eqw_dict["6.2f"])
fluxes_77 = np.array(eqw_dict["7.7f"])
fluxes_113 = np.array(eqw_dict["11.3f"])

ratio_113_77 = fluxes_113 / fluxes_77
ratio_62_77 = fluxes_62 / fluxes_77



plt.style.use('dark_background')
fig, ax = plt.subplots()
fig.set_size_inches(11, 8)
cbar=ax.scatter(ratio_62_77, ratio_113_77, c=eqw62, s=50, cmap='cool', label="Stierwalt et al. (2013)", zorder=0)
#ax.scatter(ratio_62_77[162], ratio_113_77[162], c=eqw62[162], s=100, marker="s", cmap='cool')
ax.scatter(ratio_62_77[225], ratio_113_77[225], c="#5591f2", s=300, marker="s", label="Stierwalt - IR23128", zorder=1)
#ax.scatter(ratio_62_77_new, ratio_113_77_new, c=eqw62_new, s=100, marker="*", cmap='cool')

ax.scatter(ratio_62_77_new[0], ratio_113_77_new[0], c="#a358e0", label="JWST - IR23128-S", cmap='cool', s=300, marker="v", zorder=4)
ax.scatter(ratio_62_77_new[1], ratio_113_77_new[1], c="#f72dea", label="JWST - IR23128-N", cmap='cool', s=300, marker="^", zorder=5)
#ax.scatter(ratio_62_77_new[2], ratio_113_77_new[2], c="#5591f2", label="JWST - IR15250", cmap='cool', s=300, marker="<")
ax.set_xlim(0.1, 0.5)
ax.set_ylim(0.05, 0.6)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.set_xlabel(r"L(6.2$\mu m$)/L(7.7$\mu m$" + " Complex)", fontsize=24)
ax.set_ylabel(r"L(11.3$\mu m$)/L(7.7$\mu m$" + " Complex)", fontsize=24)
plt.legend(prop={'size': 16})
colorbar = plt.colorbar(cbar, pad=0.1)
colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=24, rotation=270, labelpad=25)
colorbar.ax.yaxis.set_ticks_position("left")
colorbar.ax.tick_params(labelsize=18)
plt.savefig("./../diagnostic_plots/pah_ratios.pdf", dpi=1000)
plt.close()

plt.style.use('dark_background')
fig, ax = plt.subplots()
fig.set_size_inches(10, 10)
ax.scatter(ratio_62_77[162], ratio_113_77[162], s=100, label="Old IR15250", color="pink")
ax.scatter(ratio_62_77[225], ratio_113_77[225], s=100, label="Old IR23128", color="white")
ax.scatter(ratio_62_77_new[0], ratio_113_77_new[0], label="New IR23128-South", s=100, marker="*", color="red")
ax.scatter(ratio_62_77_new[1], ratio_113_77_new[1], label="New IR23128-North", s=100, marker="*", color="cyan")
ax.scatter(ratio_62_77_new[2], ratio_113_77_new[2], label="Old IR15250", s=100, marker="*", color="green")
ax.set_xlabel(r"L(6.2$\mu m$)/L(7.7$\mu m$" + " Complex)")
ax.set_ylabel(r"L(11.3$\mu m$)/L(7.7$\mu m$" + " Complex)")
plt.legend()
plt.savefig("./../diagnostic_plots/pah_ratios_comparison.pdf", dpi=1000)'''