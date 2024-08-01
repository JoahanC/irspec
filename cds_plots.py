from astropy.io import ascii
from astropy.table import join
from cubespec import CubeSpec

import numpy as np
import matplotlib.pyplot as plt


generate_neon_ew_plot = True
generate_pah_plot = False

st1_readme_path = "https://cdsarc.cds.unistra.fr/ftp/J/ApJS/206/1/ReadMe"
st1_file1_path = "https://cdsarc.cds.unistra.fr/ftp/J/ApJS/206/1/table1.dat"

st2_readme_path = "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/777/156/ReadMe"
st2_file1_path = "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/777/156/table1.dat"
st2_file2_path = "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/777/156/table2.dat"
st2_file3_path = "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/777/156/table3.dat"

st3_readme_path = "https://cdsarc.u-strasbg.fr/ftp/pub/J/other/ApJ/790/124/ReadMe"
st3_file1_path = "https://cdsarc.u-strasbg.fr/ftp/pub/J/other/ApJ/790/124/table1.dat"

table1 = ascii.read(st1_file1_path, readme=st1_readme_path)
table2 = ascii.read(st2_file1_path, readme=st2_readme_path)
table3 = ascii.read(st3_file1_path, readme=st3_readme_path)


if generate_neon_ew_plot:
    merged_table = join(table1, table2, keys='Name')
    merged_table_limit = merged_table[~merged_table["l_[SIV]"].mask]
    merged_table_limit = merged_table_limit[~merged_table_limit["l_[NeIII]"].mask]
    merged_table_limit = merged_table_limit[~merged_table_limit["l_[NeII]"].mask]
    merged_table_limit = merged_table_limit[~merged_table_limit["l_EW"].mask]
    merged_table_error = merged_table[~np.logical_not(merged_table["l_[SIV]"].mask)]
    merged_table_error = merged_table_error[~np.logical_not(merged_table_error["l_[NeIII]"].mask)]
    merged_table_error = merged_table_error[~np.logical_not(merged_table_error["l_[NeII]"].mask)]
    merged_table_error = merged_table_error[~np.logical_not(merged_table_error["l_EW"].mask)]

    ne3_ne2_ratio_unlog = merged_table_error["[NeIII]"] / merged_table_error["[NeII]"]
    ne3_ne2_ratio = np.log10(ne3_ne2_ratio_unlog)
    #ne3_ne2_ratio_err = ne3_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[NeIII]"] / merged_table_error["[NeIII]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    ne3_ne2_ratio_err = np.abs(np.log10(ne3_ne2_ratio_unlog) * np.log10(np.sqrt((merged_table_error["e_[NeIII]"] / merged_table_error["[NeIII]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)))
    s4_ne2_ratio_unlog = merged_table_error["[SIV]"] / merged_table_error["[NeII]"]
    s4_ne2_ratio = np.log10(s4_ne2_ratio_unlog)
    #s4_ne2_ratio_err = s4_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[SIV]"] / merged_table_error["[SIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    s4_ne2_ratio_err = np.abs(np.log10(s4_ne2_ratio_unlog) * np.log10(np.sqrt((merged_table_error["e_[SIV]"] / merged_table_error["[SIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)))
    pah_62_ew = merged_table_error["EW"]

    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    fig.set_size_inches(8, 8)
    markers, caps, bars = ax.errorbar(ne3_ne2_ratio, s4_ne2_ratio, xerr=ne3_ne2_ratio_err, yerr=s4_ne2_ratio_err, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    cvals = ax.scatter(ne3_ne2_ratio, s4_ne2_ratio, s=20, c=pah_62_ew, cmap="spring", zorder=1)
    #markers, caps, bars = ax.errorbar(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, xerr=ne3_ne2_ratio_err, yerr=s4_ne2_ratio_err, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, s=20, c=pah_62_ew, cmap="spring", zorder=1)

    colorbar = plt.colorbar(cvals)
    colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=14, rotation=270, labelpad=15)
    #[bar.set_alpha(0.5) for bar in bars]
    #[cap.set_alpha(0.5) for cap in caps]
    ax.set_xlabel("log([NeIII]/Ne[II])")
    ax.set_ylabel("log([SIV]/[NeII])")
    
    # Get new values
    
    siv_fluxes = []
    ne3_fluxes = []
    ne2_fluxes = []
    entry_idxs = [["S", "6", "AGN"], ["N", "12", "SB"], ["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"]]
    
    for entry in entry_idxs:
        spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
        gauss = spec_obj.recall_gauss()
        siv_idxs = []
        ne3_idxs = []
        ne2_idxs = []
        
        for idx, name in enumerate(gauss[3]):
            if name == "SIV_105105":
                siv_idxs.append(idx)
            if name == "NeIII_155551":
                ne3_idxs.append(idx)
            if name == "NeII_128136":
                ne2_idxs.append(idx)
        
        siv_flux = 0
        ne3_flux = 0
        ne2_flux = 0
        for idx in siv_idxs:
            siv_flux += gauss[1][idx]
        for idx in ne3_idxs:
            ne3_flux += gauss[1][idx]
        for idx in ne2_idxs:
            ne2_flux += gauss[1][idx]
        siv_fluxes.append(siv_flux)
        ne3_fluxes.append(ne3_flux)
        ne2_fluxes.append(ne2_flux)
    siv_fluxes = np.array(siv_fluxes)
    ne3_fluxes = np.array(ne3_fluxes)
    ne2_fluxes = np.array(ne2_fluxes)
    new_s4_ne2_ratio = np.log10(siv_fluxes / ne2_fluxes)
    new_ne3_ne2_ratio = np.log10(ne3_fluxes / ne2_fluxes)
    
    ax.scatter(new_ne3_ne2_ratio[0], new_s4_ne2_ratio[0], s=50, marker="*", label="IR23128S", zorder=2)
    ax.scatter(new_ne3_ne2_ratio[1], new_s4_ne2_ratio[1], s=50, marker="*", label="IR23128N", zorder=3)
    ax.scatter(new_ne3_ne2_ratio[2], new_s4_ne2_ratio[2], s=50, marker="*", label="IR23128Na", zorder=4)
    ax.scatter(new_ne3_ne2_ratio[3], new_s4_ne2_ratio[3], s=50, marker="*", label="IR23128Nb", zorder=5)
    ax.scatter(new_ne3_ne2_ratio[4], new_s4_ne2_ratio[4], s=50, marker="*", label="IR23128Nc", zorder=6)
    ax.legend()
    
    plt.show()

if generate_pah_plot:
    
    table3 = table3[~table3["EQW6.2"].mask]
    colors = ["red", "lime", "cyan", "orange"]    
    

    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    fig.set_size_inches(8, 8)
    for i in range(4):
        if i == 0:
            ntable3 = table3[table3["EQW6.2"] < 0.27]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers0, caps0, bars0 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="o", ms=5, color=colors[i], label=r"EQW$_{6.2\mu m} < 0.27$", ls="None", ecolor="white", zorder=0)
        if i == 1:
            ntable3 = table3[table3["EQW6.2"] > 0.27]
            ntable3 = ntable3[ntable3["EQW6.2"] < 0.41]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers1, caps1, bars1 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="*", ms=5, color=colors[i], label=r"0.27 < EQW$_{6.2\mu m} < 0.41$", ls="None", ecolor="white", zorder=0)
        if i == 2:
            ntable3 = table3[table3["EQW6.2"] > 0.41]
            ntable3 = ntable3[ntable3["EQW6.2"] < 0.54]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers2, caps2, bars2 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="s", ms=5, color=colors[i], label=r"0.41 < EQW$_{6.2\mu m} < 0.54$", ls="None", ecolor="white", zorder=0)
        if i == 3:
            ntable3 = table3[table3["EQW6.2"] > 0.54]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers3, caps3, bars3 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="^", ms=5, color=colors[i], label=r"EQW$_{6.2\mu m} > 0.54$", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio, s4_ne2_ratio, s=100, c=pah_62_ew, cmap="spring", zorder=1)
    #colorbar = plt.colorbar(cvals)
    #colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=14, rotation=270, labelpad=15)
    
    
    # Get new values
    
    pah11_strengths = []
    pah77_strengths = []
    pah62_strengths = []
    entry_idxs = [["S", "6", "AGN"], ["N", "12", "SB"], ["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"]]
    
    for entry in entry_idxs:
        spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
        drude = spec_obj.recall_drude()
        pah11_idxs = []
        pah77_idxs = []
        pah62_idxs = []
        
        print(drude[1])
        
        for idx, name in enumerate(drude[3]):
            if name == "9_PAH110_110200":
                pah11_idxs.append(idx)
            if name == "9_PAH112_112420":
                pah11_idxs.append(idx)
            if name == "9_PAH114_113572":
                pah11_idxs.append(idx)
            if name == "5_PAH74_74200":
                pah77_idxs.append(idx)
            if name == "5_PAH76_76000":
                pah77_idxs.append(idx)
            if name == "5_PAH78_78500":
                pah77_idxs.append(idx)
            if name == "3_PAH60_60290":
                pah62_idxs.append(idx)
            if name == "3_PAH62_62123":
                pah62_idxs.append(idx)
            if name == "3_PAH63_62700":
                pah62_idxs.append(idx)
        
        pah11_strength = 0
        pah77_strength = 0
        pah62_strength = 0
        for idx in pah11_idxs:
            pah11_strength += drude[1][idx]
        for idx in pah77_idxs:
            pah77_strength += drude[1][idx]
        for idx in pah62_idxs:
            pah62_strength += drude[1][idx]
        pah11_strengths.append(pah11_strength)
        pah77_strengths.append(pah77_strength)
        pah62_strengths.append(pah62_strength)
    

    pah11_strengths = np.array(pah11_strengths)
    

    pah77_strengths = np.array(pah77_strengths)
    pah62_strengths = np.array(pah62_strengths)
    pah_62_77_ratio = pah62_strengths / pah77_strengths 
    pah_11_77_ratio = pah11_strengths / pah77_strengths 
    print(pah_62_77_ratio[0])
    ax.scatter(pah_62_77_ratio[0], pah_11_77_ratio[0], marker="*", s=500, label="IR23128S")
    ax.scatter(pah_62_77_ratio[1], pah_11_77_ratio[1], marker="*", s=500, label="IR23128N")
    ax.scatter(pah_62_77_ratio[2], pah_11_77_ratio[2], marker="*", s=500, label="IR23128Na")
    ax.scatter(pah_62_77_ratio[3], pah_11_77_ratio[3], marker="*", s=500, label="IR23128Nb")
    ax.scatter(pah_62_77_ratio[4], pah_11_77_ratio[4], marker="*", s=500, label="IR23128Nc")
    
    
    [bar.set_alpha(0.5) for bar in bars0]
    [cap.set_alpha(0.5) for cap in caps0]
    [bar.set_alpha(0.5) for bar in bars1]
    [cap.set_alpha(0.5) for cap in caps1]
    [bar.set_alpha(0.5) for bar in bars2]
    [cap.set_alpha(0.5) for cap in caps2]
    [bar.set_alpha(0.5) for bar in bars3]
    [cap.set_alpha(0.5) for cap in caps3]
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.set_xlim(0.1, 0.55)
    #ax.set_ylim(0.05, 1.1)
    ax.legend()
    ax.set_xlabel(r"L(6.2$\mu m$)/L(7.7$\mu m$" + " Complex)")
    ax.set_ylabel(r"L(11.3$\mu m$)/L(7.7$\mu m$" + " Complex)")
    
    plt.show()




