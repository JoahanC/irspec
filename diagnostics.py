""" 
This script generates several diagnostic plots dedidated to broad 
morphological studes of IR23128 using Inami and Stierwalt values
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import join
import matplotlib.ticker as ticker


from cubespec import CubeSpec
from plotparams import PlotParams


def sulfur_neon_diagnostic(target="IR 23128-5919"):
    """ 
    Generate a radiation field hardness diagnostic using the 
    [SIV]/[NeII] ratio and the [NeIII]/[NeII] ratio.
    """
    
    # Steirwalt
    st1_readme_path = "https://cdsarc.cds.unistra.fr/ftp/J/ApJS/206/1/ReadMe"
    # Inami
    st2_readme_path = "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/777/156/ReadMe"

    # Read in both tables
    
    steirwalt_table = ascii.read(st1_file1_path, readme=st1_readme_path)
    inami_table = ascii.read(st2_file1_path, readme=st2_readme_path)
    
    ### Apply filters to remove unconstrained objects
    
    siv_ne_err = siv_ne[~np.logical_not(siv_ne["l_[SIV]"].mask)]
    siv_ne_err = siv_ne_err[~np.logical_not(siv_ne_err["l_[NeIII]"].mask)]
    siv_ne_err = siv_ne_err[~np.logical_not(siv_ne_err["l_[NeII]"].mask)]
    siv_ne_err = siv_ne_err[~np.logical_not(siv_ne_err["l_EW"].mask)]
    
    # Collect old Inami value for target 
    
    old_inami = siv_ne_err[siv_ne_err["Name"]=="IR 23128-5919"]
    
    # Calculate line ratios and conduct error propagation
    
    ne3_ne2_ratio_unlog = siv_ne_err["[NeIII]"] / siv_ne_err["[NeII]"]
    ne3_ne2_ratio = np.log10(ne3_ne2_ratio_unlog)
    ne3_ne2_ratio_err = np.abs(np.log10(ne3_ne2_ratio_unlog) * np.log10(np.sqrt((siv_ne_err["e_[NeIII]"] / siv_ne_err["[NeIII]"]) ** 2 + (siv_ne_err["e_[NeII]"] / siv_ne_err["[NeII]"]) ** 2)))
    s4_ne2_ratio_unlog = siv_ne_err["[SIV]"] / siv_ne_err["[NeII]"]
    s4_ne2_ratio = np.log10(s4_ne2_ratio_unlog)
    s4_ne2_ratio_err = np.abs(np.log10(s4_ne2_ratio_unlog) * np.log10(np.sqrt((siv_ne_err["e_[SIV]"] / siv_ne_err["[SIV]"]) ** 2 + (siv_ne_err["e_[NeII]"] / siv_ne_err["[NeII]"]) ** 2)))
    
    # Calculate Inami ratios for target 
    
    target_ne3_ne2_ratio_unlog = old_inami["[NeIII]"] / old_inami["[NeII]"]
    target_ne3_ne2_ratio = np.log10(target_ne3_ne2_ratio_unlog)
    target_s4_ne2_ratio_unlog = old_inami["[SIV]"] / old_inami["[NeII]"]
    target_s4_ne2_ratio = np.log10(target_s4_ne2_ratio_unlog)



    pass

# Instantiate plotting parameters
pltparams = PlotParams(palatte="dark", scaling="presentation")
#colors = ["tab:blue", "tab:orange", "green", "red", "purple"]
colors = ["orangered", "cyan", "magenta", "lime", "pink", "brown", "gold"]
#linestyles = pltparams.linestyles()
#markers = pltparams.markers()

# Flags for enabling which diagnostic plots to generate
generate_siv_ne_plot = True
generate_ox_iron_plot = False
generate_pah_plot = False


# Define all file paths for archival data
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


# Generates [SIV]/[NeII] vs [NeIII]/[NeII]
if generate_siv_ne_plot:
    
    # Combine U/LIRG datatable with Inami line values
    table1 = ascii.read(st1_file1_path, readme=st1_readme_path)
    table2 = ascii.read(st2_file1_path, readme=st2_readme_path)

    # Apply necessary filters to avoid data which is only upper bounded
    siv_ne = join(table1, table2, keys='Name')
    siv_ne_lim = siv_ne[~siv_ne["l_[SIV]"].mask]
    siv_ne_lim = siv_ne_lim[~siv_ne_lim["l_[NeIII]"].mask]
    siv_ne_lim = siv_ne_lim[~siv_ne_lim["l_[NeII]"].mask]
    siv_ne_lim = siv_ne_lim[~siv_ne_lim["l_EW"].mask]
    siv_ne_err = siv_ne[~np.logical_not(siv_ne["l_[SIV]"].mask)]
    siv_ne_err = siv_ne_err[~np.logical_not(siv_ne_err["l_[NeIII]"].mask)]
    siv_ne_err = siv_ne_err[~np.logical_not(siv_ne_err["l_[NeII]"].mask)]
    siv_ne_err = siv_ne_err[~np.logical_not(siv_ne_err["l_EW"].mask)]
    
    # Collect old Inami values if they exist for this object
    old_inami = siv_ne_err[siv_ne_err["Name"]=="ESO148-IG002"]

    # Calculate line ratios and conduct error propagation
    ne3_ne2_ratio_unlog = siv_ne_err["[NeIII]"] / siv_ne_err["[NeII]"]
    ne3_ne2_ratio = np.log10(ne3_ne2_ratio_unlog)
    #ne3_ne2_ratio_err = ne3_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[NeIII]"] / merged_table_error["[NeIII]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    ne3_ne2_ratio_err = np.abs(np.log10(ne3_ne2_ratio_unlog) * np.log10(np.sqrt((siv_ne_err["e_[NeIII]"] / siv_ne_err["[NeIII]"]) ** 2 + (siv_ne_err["e_[NeII]"] / siv_ne_err["[NeII]"]) ** 2)))
    s4_ne2_ratio_unlog = siv_ne_err["[SIV]"] / siv_ne_err["[NeII]"]
    s4_ne2_ratio = np.log10(s4_ne2_ratio_unlog)
    #s4_ne2_ratio_err = s4_ne2_ratio_unlog * np.sqrt((merged_table_error["e_[SIV]"] / merged_table_error["[SIV]"]) ** 2 + (merged_table_error["e_[NeII]"] / merged_table_error["[NeII]"]) ** 2)
    s4_ne2_ratio_err = np.abs(np.log10(s4_ne2_ratio_unlog) * np.log10(np.sqrt((siv_ne_err["e_[SIV]"] / siv_ne_err["[SIV]"]) ** 2 + (siv_ne_err["e_[NeII]"] / siv_ne_err["[NeII]"]) ** 2)))
    pah_62_ew = siv_ne_err["EW"]
    
    # Calculate ratio for old Inami values of IR23128
    eso_ne3_ne2_ratio_unlog = old_inami["[NeIII]"] / old_inami["[NeII]"]
    eso_ne3_ne2_ratio = np.log10(eso_ne3_ne2_ratio_unlog)
    eso_s4_ne2_ratio_unlog = old_inami["[SIV]"] / old_inami["[NeII]"]
    eso_s4_ne2_ratio = np.log10(eso_s4_ne2_ratio_unlog)

    

    #markers, caps, bars = ax.errorbar(ne3_ne2_ratio, s4_ne2_ratio, xerr=ne3_ne2_ratio_err * 0.05, yerr=s4_ne2_ratio_err * 0.05, alpha=0.3, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio, s4_ne2_ratio, s=20, c=pah_62_ew, alpha=0.8, cmap="spring", zorder=1, label="Inami+2013")
    
    #markers, caps, bars = ax.errorbar(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, xerr=ne3_ne2_ratio_err, yerr=s4_ne2_ratio_err, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, s=20, c=pah_62_ew, cmap="spring", zorder=1)

    #colorbar = plt.colorbar(cvals)
    #colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=18, rotation=270, labelpad=15)
    #colorbar.ax.tick_params(labelsize=16)

    #[bar.set_alpha(0.5) for bar in bars]
    #[cap.set_alpha(0.5) for cap in caps]
    
    
    
    
    # Get new values
    
    siv_fluxes = []
    ne3_fluxes = []
    ne2_fluxes = []
    entry_idxs = [["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"], ["N", "12", "SB"], ["S", "6", "AGN"]]
    

    for entry in entry_idxs:
        spec_obj = CubeSpec("./../", "param_files", f"IR23128-{entry[0]}_{entry[1]}_single_param.txt", f"input_data/IR23128-{entry[0]}/", redshift=0.044601, fit_dirname=f"IR23128{entry[0]}_{entry[2]}{entry[1]}", mode=entry[2])
        filepath = spec_obj.return_lines()
        pahfile = spec_obj.return_pah()
        data = spec_obj.recall_data()
        with open(filepath, 'r') as filey:
            data = filey.readlines()
        siv_flux = 0
        ne3_flux = 0
        ne2_flux = 0
        for line in data:
            if "SIV_105105" in line:
                siv_flux += float(line.split()[4])
            if "NeIII_155551" in line:
                ne3_flux += float(line.split()[4])
            if "NeII_128136" in line:
                pass
                #ne2_flux += float(line.split()[4])
        
        #siv_fluxes.append(siv_flux)
        siv_fluxes.append(spec_obj.local_fit_flux([10.45, 10.575], 10.507, "[SIV]", 2, 2).value)
        ne2_fluxes.append(spec_obj.local_fit_flux([12.72, 12.84], 12.813, "[NeII]", 1, 2).value)
        
        spec_obj.local_fit_flux([12.72, 12.84], 12.813, "[NeII]", 1, 2)
        spec_obj.local_fit_flux([15.39, 15.75], 15.556, "[NeIII]", 1, 3)
        spec_obj.local_fit_flux([10.45, 10.575], 10.507, "[SIV]", 1, 2)
        
        if entry[1] in ["0", "1", "2"]:
            ne3_fluxes.append(spec_obj.local_fit_flux([15.39, 15.75], 15.556, "[NeIII]", 1, 2).value)
        if entry[1] in ["6", "12"]:
            ne3_fluxes.append(spec_obj.local_fit_flux([15.39, 15.75], 15.556, "[NeIII]", 1, 3).value)
        #ne3_fluxes.append(ne3_flux)
        #ne2_fluxes.append(spec_obj.local_fit_flux([12.5, 13.1], 12.813, "[NeII]", 1, 2).value)
        
        #ne2_fluxes.append(ne2_flux)
    
    pltparams = PlotParams(palatte="dark", scaling="presentation")

    
    #plt.style.use('default')
    fig, ax = plt.subplots()
    fig.set_size_inches(16, 10)
    
    #markers, caps, bars = ax.errorbar(s4_ne2_ratio, ne3_ne2_ratio, xerr=s4_ne2_ratio_err * 0.1, yerr=ne3_ne2_ratio_err * 0.05, alpha=0.3, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    cvals = ax.scatter(s4_ne2_ratio, ne3_ne2_ratio, s=100, c="white", alpha=0.4, zorder=1, label="Inami+2013")
    ax.set_xlabel("[NeIII]/[NeII]")
    ax.set_ylabel("[SIV]/[NeII]")
    
    siv_fluxes = np.array(siv_fluxes)
    ne3_fluxes = np.array(ne3_fluxes)
    ne2_fluxes = np.array(ne2_fluxes)
    new_s4_ne2_ratio = np.log10(siv_fluxes / ne2_fluxes)
    new_ne3_ne2_ratio = np.log10(ne3_fluxes / ne2_fluxes)
    print(new_s4_ne2_ratio)
    
    #colors = ["tab:blue", "tab:orange", "green", "gold", "purple"]
    ["orangered", "cyan", "magenta", "lime", "pink", "brown", "gold"]
    colors = pltparams.dark_colors()


    # Plot all of the new data points

    """ax.scatter(new_ne3_ne2_ratio[0], new_s4_ne2_ratio[0], s=800, marker="v", label=r"S", zorder=2)
    ax.scatter(new_ne3_ne2_ratio[1], new_s4_ne2_ratio[1], s=800, marker="^", label=r"N", zorder=3)
    ax.scatter(new_ne3_ne2_ratio[2], new_s4_ne2_ratio[2], s=800, marker="*", label=r"Na", zorder=4)
    ax.scatter(new_ne3_ne2_ratio[3], new_s4_ne2_ratio[3], s=800, marker="*", label=r"Nb", zorder=5)
    ax.scatter(new_ne3_ne2_ratio[4], new_s4_ne2_ratio[4], s=800, marker="*", label=r"Nc", zorder=6)"""
    
    # Plot all of the new data points
    ax.scatter(new_s4_ne2_ratio[0], new_ne3_ne2_ratio[0], s=1000, marker="*", color=colors[0], label=r"a", edgecolors="black", zorder=6)
    ax.scatter(new_s4_ne2_ratio[1], new_ne3_ne2_ratio[1], s=1000, marker="*", color=colors[1], label=r"b", edgecolors="black", zorder=5)
    ax.scatter(new_s4_ne2_ratio[2], new_ne3_ne2_ratio[2], s=1000, marker="*", color=colors[2], label=r"c", edgecolors="black", zorder=4)
    ax.scatter(new_s4_ne2_ratio[3], new_ne3_ne2_ratio[3], s=1000, marker="^", color=colors[3], label=r"N", edgecolors="black", zorder=3)
    ax.scatter(new_s4_ne2_ratio[4], new_ne3_ne2_ratio[4], s=1000, marker="v", color=colors[4], label=r"S", edgecolors="black", zorder=2)
    ax.scatter(eso_s4_ne2_ratio[0], eso_ne3_ne2_ratio[0], s=1000, color=colors[5], label="IR 23128-5919", zorder=1)

    # Plot old values
    
    
    """# Plot old values
    ax.scatter(eso_ne3_ne2_ratio[0], eso_s4_ne2_ratio[0], s=200, color=colors[3], zorder=8, label="IR23128-Inami")"""

    #ax.legend()
    ax.set_xlabel("log([SIV]/[NeII])", fontsize=32)
    ax.set_ylabel("log([NeIII]/[NeII])", fontsize=32)
    ax.tick_params(axis='x', labelsize=28)
    ax.tick_params(axis='y', labelsize=28)
    ax.tick_params(direction='in', which='both', length=6, width=1, top=True)
    #plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax.get_yticklabels(), visible=False)
    plt.legend(prop={'size': 24})
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))

    plt.savefig("./../diagnostic_plots/lines/siv_ne2_vs_ne3_ne1.pdf", dpi=1200, bbox_inches="tight")


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


    #plt.style.use('dark_background')
    fig, ax = plt.subplots(facecolor="white")
    fig.set_size_inches(16, 10)

    #fig.set_size_inches(8, 8)
    #markers, caps, bars = ax.errorbar(o4_ne2_ratio, fe2_o4_ratio, xerr=o4_ne2_ratio_err * 0.05, yerr=fe2_o4_ratio_err * 0.05, alpha=0.8, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    
    #markers, caps, bars = ax.errorbar(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, xerr=ne3_ne2_ratio_err, yerr=s4_ne2_ratio_err, fmt="o", ms=0, color="white", ls="None", ecolor="white", zorder=0)
    #cvals = ax.scatter(ne3_ne2_ratio_unlog, s4_ne2_ratio_unlog, s=20, c=pah_62_ew, cmap="spring", zorder=1)

    #colorbar = plt.colorbar(cvals)
    #colorbar.set_label(r"6.2$\mu m$ PAH EW", fontsize=18, rotation=270, labelpad=15)
    #colorbar.ax.tick_params(labelsize=16)

    #[bar.set_alpha(0.5) for bar in bars]
    #[cap.set_alpha(0.5) for cap in caps]
    ax.set_xlabel("[OIV]/[NeII]", fontsize=32)
    ax.set_ylabel("[FeII]/[OIV]", fontsize=32)
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
                oiv_flux += float(line.split()[4])
            if "FeII_259883" in line:
                fe2_flux += float(line.split()[4])
            #if "NeII_128136" in line:
            #    ne2_flux += float(line.split()[4])
        data = spec_obj.recall_data()
        
        """
        gauss = spec_obj.return_gauss()
        for idx, name in enumerate(gauss[3]):
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
        #ne2_fluxes.append(ne2_flux)
        ne2_fluxes.append(spec_obj.local_fit_flux([12.5, 13.1], 12.813, "[NeII]", 1, 2).value)
    
    pltparams = PlotParams(palatte="dark", scaling="presentation")

    fig, ax = plt.subplots(facecolor="white")
    fig.set_size_inches(10, 14)
    
    ax.set_xlabel("[OIV]/[NeII]", fontsize=32)
    ax.set_ylabel("[FeII]/[OIV]", fontsize=32)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=18)
    
    cvals = ax.scatter(o4_ne2_ratio, fe2_o4_ratio, color="black", s=200, alpha=0.8, zorder=1, label="Inami+2013")
    
    oiv_fluxes = np.array(oiv_fluxes)
    fe2_fluxes = np.array(fe2_fluxes)
    ne2_fluxes = np.array(ne2_fluxes)
    new_o4_ne2_ratio = np.log10(oiv_fluxes / ne2_fluxes)
    new_fe2_o4_ratio = np.log10(fe2_fluxes / oiv_fluxes)
    #print(new_s4_ne2_ratio)
    
    #ax.scatter(new_o4_ne2_ratio[0], new_fe2_o4_ratio[0], s=500, marker="v", color=colors[0], label="S", zorder=2)
    #ax.scatter(new_o4_ne2_ratio[1], new_fe2_o4_ratio[1], s=500, marker="^", color=colors[1], label="N", zorder=3)
    ax.scatter(new_o4_ne2_ratio[0], new_fe2_o4_ratio[0], s=500, marker="v", label="S", zorder=2)
    ax.scatter(new_o4_ne2_ratio[1], new_fe2_o4_ratio[1], s=500, marker="^", label="N", zorder=3)
    #ax.scatter(new_o4_ne2_ratio[2], new_fe2_o4_ratio[2], s=200, marker=markers[2], color=colors[2], label="Na", zorder=4)
    #ax.scatter(new_o4_ne2_ratio[3], new_fe2_o4_ratio[3], s=200, marker=markers[2], color=colors[3], label="Nb", zorder=5)
    #ax.scatter(new_o4_ne2_ratio[4], new_fe2_o4_ratio[4], s=200, marker=markers[2], color=colors[4], label="Nc", zorder=6)
    
    #ax.scatter(eso_o4_ne2_ratio[0], eso_fe2_o4_ratio[0], s=500, color=colors[5], label="IR23128-Inami", zorder=6)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.tick_params(direction='in', which='both', length=6, width=1, top=True)
    #plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax.get_yticklabels(), visible=False)
    plt.legend(prop={'size': 24})
    ax.tick_params(direction='in', length=6, width=2, grid_alpha=0.5)

    #plt.show()
    plt.savefig("./../diagnostic_plots/oxiron_pos.pdf", dpi=1200, bbox_inches="tight")



if generate_pah_plot:
    
    table3 = table3[~table3["EQW6.2"].mask]
    #colors = ["red", "lime", "cyan", "orange"]    
    

    #plt.style.use('dark_background')
    #fig, ax = plt.subplots()

    #fig.set_size_inches(12, 8)
    
    pah11_pah7_ratio = table3["11.3PAH"] / table3["7.7PAH"]
    pah11_pah7_ratio_err = pah11_pah7_ratio * np.sqrt((table3["e_11.3PAH"] / table3["11.3PAH"]) ** 2 + (table3["e_7.7PAH"] / table3["7.7PAH"]) ** 2)
    pah6_pah7_ratio = table3["6.2PAH"] / table3["7.7PAH"]
    pah6_pah7_ratio_err = pah6_pah7_ratio * np.sqrt((table3["e_6.2PAH"] / table3["6.2PAH"]) ** 2 + (table3["e_6.2PAH"] / table3["6.2PAH"]) ** 2)
    
    
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
    
    
    # Get new values
    
    pah11_strengths = []
    pah77_strengths = []
    pah62_strengths = []
    entry_idxs = [["N", "0", "SB"], ["N", "1", "SB"], ["N", "2", "SB"], ["N", "12", "SB"], ["S", "6", "AGN"]]
    
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
        #print(entry, pah62_strength)
        
        pah11_strengths.append(pah11_strength)
        pah77_strengths.append(pah77_strength)
        pah62_strengths.append(pah62_strength)
    
    
    pltparams = PlotParams(palatte="dark", scaling="presentation")
    colors = pltparams.dark_colors()
    fig, ax = plt.subplots()
    fig.set_size_inches(16, 10)
    
    colors = ["orangered", "cyan", "magenta", "lime", "pink", "brown", "gold"]
    colors = pltparams.dark_colors()


    ax.scatter(pah6_pah7_ratio, pah11_pah7_ratio, label="Stierwalt+2014", s=100, color="white", alpha=0.4, zorder=1)
    markers0, caps0, bars0 = ax.errorbar(pah6_pah7_ratio, pah11_pah7_ratio, xerr=pah6_pah7_ratio_err, yerr=pah11_pah7_ratio_err, fmt="o", ms=2, color="white", ls="None", ecolor="white", alpha=0.3, zorder=0)


    pah11_strengths = np.array(pah11_strengths)
    

    pah77_strengths = np.array(pah77_strengths)
    pah62_strengths = np.array(pah62_strengths)
    pah_62_77_ratio = pah62_strengths / pah77_strengths 
    pah_11_77_ratio = pah11_strengths / pah77_strengths 
    print(colors)

    #ax.scatter(pah_62_77_ratio[0], pah_11_77_ratio[0], s=800, marker="*", color=colors[0], label="S")
    #ax.scatter(pah_62_77_ratio[1], pah_11_77_ratio[1], s=800, marker="*", color=colors[1], label="N")
    #ax.scatter(pah_62_77_ratio[2], pah_11_77_ratio[2], s=800, marker="*", color=colors[2], label="Na")
    #ax.scatter(pah_62_77_ratio[3], pah_11_77_ratio[3], s=800, marker="^", color=colors[3], label="Nb")
    #ax.scatter(pah_62_77_ratio[4], pah_11_77_ratio[4], s=800, marker="v", color=colors[4], label="Nc")
    ax.scatter(pah_62_77_ratio[0], pah_11_77_ratio[0], s=1000, marker="*", color=colors[0], label="a", edgecolor="black")
    ax.scatter(pah_62_77_ratio[1], pah_11_77_ratio[1], s=1000, marker="*", color=colors[1], label="b", edgecolor="black")
    ax.scatter(pah_62_77_ratio[2], pah_11_77_ratio[2], s=1000, marker="*", color=colors[2], label="c", edgecolor="black")
    ax.scatter(pah_62_77_ratio[3], pah_11_77_ratio[3], s=1000, marker="v", color=colors[3], label="S", edgecolor="black")
    ax.scatter(pah_62_77_ratio[4], pah_11_77_ratio[4], s=1000, marker="^", color=colors[4], label="N", edgecolor="black")
    #ax.scatter(pah_62_77_ratio[4], pah_11_77_ratio[4], s=200, marker=markers[2], color=colors[4], label="Nc")
    #ax.scatter(old_6_7_ratio, old_11_7_ratio, s=200, marker=markers[3], color=colors[5], label="IR23128-Stierwalt")
    
    

    
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
    ax.tick_params(direction='in', which='both', length=6, width=0, top=True)
    #plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax.get_yticklabels(), visible=False)
    #plt.legend(prop={'size': 24})
    
    ax.legend(prop={'size': 24})
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))

    ax.set_xlabel(r"L(6.2$\mu m$)/L(7.7$\mu m$" + " Complex)", fontsize=32)
    ax.set_ylabel(r"L(11.3$\mu m$)/L(7.7$\mu m$" + " Complex)", fontsize=32)
    ax.tick_params(axis='both', which='major', labelsize=28)
    #ax.tick_params(axis='both', which='minor', labelsize=18)
    ax.tick_params(direction='in', length=6, width=1, grid_alpha=0.5)
    
    #plt.show()
    plt.savefig("./../diagnostic_plots/pah_ionization.pdf", dpi=1200, bbox_inches="tight")




