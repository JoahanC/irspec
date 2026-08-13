from astropy.io import ascii
from astropy.table import join
from cubespec import CubeSpec

import numpy as np
import matplotlib.pyplot as plt


generate_neon_ew_plot = False
generate_ox_iron_plot = False
generate_pah_plot = True

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
    
    old_inami = merged_table[merged_table["Name"]=="ESO148-IG002"]

    ne3_ne2_ratio_unlog = merged_table_error["[NeIII]"] / merged_table_error["[NeII]"]
    ne3_ne2_ratio = np.log10(ne3_ne2_ratio_unlog)
    #ne3_ne2_ratio_err = ne3_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[NeIII]"] / merged_table_error["[NeIII]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    ne3_ne2_ratio_err = np.abs(np.log10(ne3_ne2_ratio_unlog) * np.log10(np.sqrt((merged_table_error["e_[NeIII]"] / merged_table_error["[NeIII]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)))
    s4_ne2_ratio_unlog = merged_table_error["[SIV]"] / merged_table_error["[NeII]"]
    s4_ne2_ratio = np.log10(s4_ne2_ratio_unlog)
    #s4_ne2_ratio_err = s4_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[SIV]"] / merged_table_error["[SIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    s4_ne2_ratio_err = np.abs(np.log10(s4_ne2_ratio_unlog) * np.log10(np.sqrt((merged_table_error["e_[SIV]"] / merged_table_error["[SIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)))
    pah_62_ew = merged_table_error["EW"]
    
    eso_ne3_ne2_ratio_unlog = old_inami["[NeIII]"] / old_inami["[NeII]"]
    eso_ne3_ne2_ratio = np.log10(eso_ne3_ne2_ratio_unlog)
    eso_s4_ne2_ratio_unlog = old_inami["[SIV]"] / old_inami["[NeII]"]
    eso_s4_ne2_ratio = np.log10(eso_s4_ne2_ratio_unlog)

    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    fig.set_size_inches(8, 8)
    markers, caps, bars = ax.errorbar(ne3_ne2_ratio, s4_ne2_ratio, xerr=ne3_ne2_ratio_err * 0.05, yerr=s4_ne2_ratio_err * 0.05, alpha=0.3, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    cvals = ax.scatter(ne3_ne2_ratio, s4_ne2_ratio, s=20, c=pah_62_ew, alpha=0.8, cmap="spring", zorder=1, label="Inami+2013")
    #markers, caps, bars = ax.errorbar(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, xerr=ne3_ne2_ratio_err, yerr=s4_ne2_ratio_err, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, s=20, c=pah_62_ew, cmap="spring", zorder=1)

    colorbar = plt.colorbar(cvals)
    colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=18, rotation=270, labelpad=15)
    colorbar.ax.tick_params(labelsize=16)

    #[bar.set_alpha(0.5) for bar in bars]
    #[cap.set_alpha(0.5) for cap in caps]
    ax.set_xlabel("log([NeIII]/[NeII])", fontsize=18)
    ax.set_ylabel("log([SIV]/[NeII])", fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=18)
    
    ax.scatter(eso_ne3_ne2_ratio[0], eso_s4_ne2_ratio[0], s=200, marker="^", zorder=8, label="IR23128-Inami")
    
    # Get new values
    
    siv_fluxes = []
    ne3_fluxes = []
    ne2_fluxes = []
    entry_idxs = [["S", "6", "AGN"], ["N", "12", "SB"], ["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"]]
    
    for entry in entry_idxs:
        spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
        filepath = spec_obj.return_lines()
        pahfile = spec_obj.return_pah()
        with open(filepath, 'r') as filey:
            data = filey.readlines()
        siv_flux = 0
        ne3_flux = 0
        ne2_flux = 0
        for line in data:
            if "SIV_105105" in line:
                siv_flux += float(line.split()[2])
            if "NeIII_155551" in line:
                ne3_flux += float(line.split()[2])
            if "NeII_128136" in line:
                ne2_flux += float(line.split()[2])
        
        siv_fluxes.append(siv_flux)
        ne3_fluxes.append(ne3_flux)
        ne2_fluxes.append(ne2_flux)
    siv_fluxes = np.array(siv_fluxes)
    ne3_fluxes = np.array(ne3_fluxes)
    ne2_fluxes = np.array(ne2_fluxes)
    new_s4_ne2_ratio = np.log10(siv_fluxes / ne2_fluxes)
    new_ne3_ne2_ratio = np.log10(ne3_fluxes / ne2_fluxes)
    print(new_s4_ne2_ratio)
    
    
    ax.scatter(new_ne3_ne2_ratio[0], new_s4_ne2_ratio[0], s=200, marker="*", label=r"S", zorder=2)
    ax.scatter(new_ne3_ne2_ratio[1], new_s4_ne2_ratio[1], s=200, marker="*", label=r"N", zorder=3)
    ax.scatter(new_ne3_ne2_ratio[2], new_s4_ne2_ratio[2], s=200, marker="*", label=r"Na", zorder=4)
    ax.scatter(new_ne3_ne2_ratio[3], new_s4_ne2_ratio[3], s=200, marker="*", label=r"Nb", zorder=5)
    ax.scatter(new_ne3_ne2_ratio[4], new_s4_ne2_ratio[4], s=200, marker="*", label=r"Nc", zorder=6)

    ax.legend(prop={'size': 14})

    #plt.show()
    plt.savefig("./emissionline.pdf", dpi=1000, bbox_inches="tight")


if generate_ox_iron_plot:
    table4 = ascii.read(st2_file2_path, readme=st2_readme_path)
    merged_table = join(table1, table2, keys='Name')
    merged_table = join(merged_table, table4, keys='Name')
    merged_table_limit = merged_table[~merged_table["l_[OIV]"].mask]
    merged_table_limit = merged_table_limit[~merged_table_limit["l_[FeII]"].mask]
    merged_table_limit = merged_table_limit[~merged_table_limit["l_[NeII]"].mask]
    merged_table_error = merged_table[~np.logical_not(merged_table["l_[OIV]"].mask)]
    merged_table_error = merged_table_error[~np.logical_not(merged_table_error["l_[FeII]"].mask)]
    merged_table_error = merged_table_error[~np.logical_not(merged_table_error["l_[NeII]"].mask)]
    
    old_inami = merged_table[merged_table["Name"]=="ESO148-IG002"]

    fe2_o4_ratio_unlog = merged_table_error["[FeII]"] / merged_table_error["[OIV]"]
    fe2_o4_ratio = np.log10(fe2_o4_ratio_unlog)
    #ne3_ne2_ratio_err = ne3_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[NeIII]"] / merged_table_error["[NeIII]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    fe2_o4_ratio_err = np.abs(np.log10(fe2_o4_ratio_unlog) * np.log10(np.sqrt((merged_table_error["e_[FeII]"] / merged_table_error["[FeII]"]) ** 2 + (merged_table_error["e_[OIV]"] / merged_table_error["[OIV]"]) ** 2)))
    o4_ne2_ratio_unlog = merged_table_error["[OIV]"] / merged_table_error["[NeII]"]
    o4_ne2_ratio = np.log10(o4_ne2_ratio_unlog)
    #s4_ne2_ratio_err = s4_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[SIV]"] / merged_table_error["[SIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    o4_ne2_ratio_err = np.abs(np.log10(o4_ne2_ratio_unlog) * np.log10(np.sqrt((merged_table_error["e_[OIV]"] / merged_table_error["[OIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)))

    eso_fe2_o4_ratio_unlog = old_inami["[FeII]"] / old_inami["[OIV]"]
    eso_fe2_o4_ratio = np.log10(eso_fe2_o4_ratio_unlog)
    eso_o4_ne2_ratio_unlog = old_inami["[OIV]"] / old_inami["[NeII]"]
    eso_o4_ne2_ratio = np.log10(eso_o4_ne2_ratio_unlog)


    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    fig.set_size_inches(8, 8)
    markers, caps, bars = ax.errorbar(o4_ne2_ratio, fe2_o4_ratio, xerr=o4_ne2_ratio_err * 0.05, yerr=fe2_o4_ratio_err * 0.05, alpha=0.8, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    cvals = ax.scatter(o4_ne2_ratio, fe2_o4_ratio, s=20, alpha=0.8, zorder=1, color="cyan",label="Inami+2013")
    #markers, caps, bars = ax.errorbar(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, xerr=ne3_ne2_ratio_err, yerr=s4_ne2_ratio_err, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, s=20, c=pah_62_ew, cmap="spring", zorder=1)

    #colorbar = plt.colorbar(cvals)
    #colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=18, rotation=270, labelpad=15)
    #colorbar.ax.tick_params(labelsize=16)

    #[bar.set_alpha(0.5) for bar in bars]
    #[cap.set_alpha(0.5) for cap in caps]
    ax.set_xlabel("log([OIV]/[FeII])", fontsize=18)
    ax.set_ylabel("log([FeII]/[OIV])", fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=18)
    
    #ax.scatter(eso_o4_ne2_ratio[0], eso_fe2_o4_ratio[0], s=200, marker="^", zorder=8, label="IR23128-Inami")
    
    # Get new values
    
    oiv_fluxes = []
    fe2_fluxes = []
    ne2_fluxes = []
    entry_idxs = [["S", "6", "AGN"], ["N", "12", "SB"]]
    
    for entry in entry_idxs:
        spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
        filepath = spec_obj.return_lines()
        
        with open(filepath, 'r') as filey:
            data = filey.readlines()
        
        oiv_flux = 0
        fe2_flux = 0
        ne2_flux = 0
        
        for line in data:
            if "OIV_258903" in line:
                oiv_flux += float(line.split()[2])
            if "FeII_259883" in line:
                fe2_flux += float(line.split()[2])
            if "NeII_128136" in line:
                ne2_flux += float(line.split()[2])
        
        """for idx, name in enumerate(gauss[3]):
            if name == "OIV_258903":
                oiv_idxs.append(idx)
            if name == "FeII_259883":
                fe2_idxs.append(idx)
            if name == "NeII_128136":
                ne2_idxs.append(idx)
        
        for idx in oiv_idxs:
            oiv_flux += gauss[2][idx]
        for idx in fe2_idxs:
            fe2_flux += gauss[2][idx]
        for idx in ne2_idxs:
            ne2_flux += gauss[2][idx]"""
        
        oiv_fluxes.append(oiv_flux)
        fe2_fluxes.append(fe2_flux)
        ne2_fluxes.append(ne2_flux)
    oiv_fluxes = np.array(oiv_fluxes)
    fe2_fluxes = np.array(fe2_fluxes)
    ne2_fluxes = np.array(ne2_fluxes)
    new_o4_ne2_ratio = np.log10(oiv_fluxes / ne2_fluxes)
    new_fe2_o4_ratio = np.log10(fe2_fluxes / oiv_fluxes)
    #print(new_s4_ne2_ratio)
    
    ax.scatter(new_o4_ne2_ratio[0], new_fe2_o4_ratio[0], s=200, marker="*", label=r"South", zorder=2, color="lime")
    ax.scatter(new_o4_ne2_ratio[1], new_fe2_o4_ratio[1], s=200, marker="*", label=r"North", zorder=3, color="orangered")
    #ax.scatter(new_ne3_ne2_ratio[2], new_s4_ne2_ratio[2], s=200, marker="*", label=r"Na, PAH$_{6.2\mu m} = 2.3$", zorder=4)
    #ax.scatter(new_ne3_ne2_ratio[3], new_s4_ne2_ratio[3], s=200, marker="*", label=r"Nb, PAH$_{6.2\mu m} = 1.2$", zorder=5)
    #ax.scatter(new_ne3_ne2_ratio[4], new_s4_ne2_ratio[4], s=200, marker="*", label=r"Nc, PAH$_{6.2\mu m} = 1.9$", zorder=6)
    ax.legend(prop={'size': 14})

    #plt.show()
    plt.savefig("./oxiron.pdf", dpi=1000, bbox_inches="tight")



if generate_pah_plot:
    
    table3 = table3[~table3["EQW6.2"].mask]
    colors = ["red", "lime", "cyan", "orange"]    
    

    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    fig.set_size_inches(12, 8)
    
    pah11_pah7_ratio = table3["11.3PAH"] / table3["7.7PAH"]
    pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((table3["e_11.3PAH"] / table3["11.3PAH"]) ** 2 + (table3["e_7.7PAH"] / table3["7.7PAH"]) ** 2)
    pah6_pah7_ratio = table3["6.2PAH"] / table3["7.7PAH"]
    pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((table3["e_6.2PAH"] / table3["6.2PAH"]) ** 2 + (table3["e_6.2PAH"] / table3["6.2PAH"]) ** 2)
    ax.scatter(pah6_pah7_ratio, pah11_pah7_ratio, label="Stierwalt+2014", s=15, c="red", zorder=1)
    markers0, caps0, bars0 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="o", ms=2, color="white", ls="None", ecolor="white", alpha=0.5, zorder=0)
    
    """for i in range(4):
        if i == 0:
            ntable3 = table3[table3["EQW6.2"] < 0.27]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers0, caps0, bars0 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="o", ms=5, color=colors[i], label=r"EQW$_{6.2\mu m} < 0.27$", ls="None", ecolor="white", alpha=0.5, zorder=0)
        if i == 1:
            ntable3 = table3[table3["EQW6.2"] > 0.27]
            ntable3 = ntable3[ntable3["EQW6.2"] < 0.41]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers1, caps1, bars1 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="*", ms=5, color=colors[i], label=r"0.27 < EQW$_{6.2\mu m} < 0.41$", ls="None", ecolor="white", alpha=1, zorder=0)
        if i == 2:
            ntable3 = table3[table3["EQW6.2"] > 0.41]
            ntable3 = ntable3[ntable3["EQW6.2"] < 0.54]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers2, caps2, bars2 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="s", ms=5, color=colors[i], label=r"0.41 < EQW$_{6.2\mu m} < 0.54$", ls="None", ecolor="white", alpha=1, zorder=0)
        if i == 3:
            ntable3 = table3[table3["EQW6.2"] > 0.54]
            pah11_pah7_ratio = ntable3["11.3PAH"] / ntable3["7.7PAH"]
            pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((ntable3["e_11.3PAH"] / ntable3["11.3PAH"]) ** 2 + (ntable3["e_7.7PAH"] / ntable3["7.7PAH"]) ** 2)
            pah6_pah7_ratio = ntable3["6.2PAH"] / ntable3["7.7PAH"]
            pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2 + (ntable3["e_6.2PAH"] / ntable3["6.2PAH"]) ** 2)
            markers3, caps3, bars3 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="^", ms=5, color=colors[i], label=r"EQW$_{6.2\mu m} > 0.54$", ls="None", ecolor="white", alpha=1, zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio, s4_ne2_ratio, s=100, c=pah_62_ew, cmap="spring", zorder=1)"""
    #colorbar = plt.colorbar(cvals)
    #colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=14, rotation=270, labelpad=15)
    old_steir = table3[table3["Name"]=="ESO148-IG002"]
    old_11_7_ratio = old_steir["11.3PAH"] / old_steir["7.7PAH"]
    old_6_7_ratio = old_steir["6.2PAH"] / old_steir["7.7PAH"]
    ax.scatter(old_6_7_ratio, old_11_7_ratio, label="Stierwalt-IR23128", s=500, marker="^")
    
    # Get new values
    
    pah11_strengths = []
    pah77_strengths = []
    pah62_strengths = []
    entry_idxs = [["S", "6", "AGN"], ["N", "12", "SB"], ["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"]]
    
    for entry in entry_idxs:
        spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
        drude = spec_obj.return_pah()
        
        with open(drude, 'r') as filey:
            data = filey.readlines()
        
        pah11_strength = 0
        pah77_strength = 0
        pah62_strength = 0
        
        for line in data:
            if "PAH62" in line:
                pah62_strength += float(line.split()[1])
            if "PAH77_C" in line:
                pah77_strength += float(line.split()[1])
            if "PAH113_C" in line:
                pah11_strength += float(line.split()[1])
        
        pah11_strengths.append(pah11_strength)
        pah77_strengths.append(pah77_strength)
        pah62_strengths.append(pah62_strength)
    

    pah11_strengths = np.array(pah11_strengths)
    

    pah77_strengths = np.array(pah77_strengths)
    pah62_strengths = np.array(pah62_strengths)
    pah_62_77_ratio = pah62_strengths / pah77_strengths 
    pah_11_77_ratio = pah11_strengths / pah77_strengths 
    #print(pah_62_77_ratio[0])
    ax.scatter(pah_62_77_ratio[0], pah_11_77_ratio[0], marker="*", s=500, label="S")
    ax.scatter(pah_62_77_ratio[1], pah_11_77_ratio[1], marker="*", s=500, label="N")
    ax.scatter(pah_62_77_ratio[2], pah_11_77_ratio[2], marker="*", s=500, label="Na")
    ax.scatter(pah_62_77_ratio[3], pah_11_77_ratio[3], marker="*", s=500, label="Nb")
    ax.scatter(pah_62_77_ratio[4], pah_11_77_ratio[4], marker="*", s=500, label="Nc")
    
    ax.set_xlim(0.1, 0.5)
    ax.set_ylim(0, 0.7)
    [bar.set_alpha(0.5) for bar in bars0]
    [cap.set_alpha(0.5) for cap in caps0]
    """[bar.set_alpha(0.5) for bar in bars1]
    [cap.set_alpha(0.5) for cap in caps1]
    [bar.set_alpha(0.5) for bar in bars2]
    [cap.set_alpha(0.5) for cap in caps2]
    [bar.set_alpha(0.5) for bar in bars3]
    [cap.set_alpha(0.5) for cap in caps3]"""
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    #ax.set_xlim(0.1, 0.55)
    #ax.set_ylim(0.05, 0.6)
    ax.legend()
    ax.set_xlabel(r"L(6.2$\mu m$)/L(7.7$\mu m$" + " Complex)", fontsize=18)
    ax.set_ylabel(r"L(11.3$\mu m$)/L(7.7$\mu m$" + " Complex)", fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=18)
    
    #plt.show()
    plt.savefig("./phaplot.pdf", dpi=1000)




