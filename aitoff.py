import numpy as np 
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


filepath = "./../archival_data/goals_table.tbl"
datatable = Table.read(filepath, format="ipac")

obj_ras = datatable["RA"]
obj_decs = datatable["Dec"]
obj_redshifts = datatable["Redshift (z)"]
print(obj_redshifts)

fixed_names = []
for name in datatable["Object Name"]:
    new_name = name.replace('<a HREF="#" onclick="getObjectSearchInfo(jQuery, 2, this); return false;">', '')
    fixed_names.append(new_name.replace("</a>", ""))


ulirgs = ["WISEA J052101.39-252145.2", "WISEA J072737.58-025454.1", "IRAS 08572+3915", "WISEA J090412.70-362701.5", "WISEA J091338.84-101920.0", "UGC 05101",
          "WISEA J105918.13+243234.5", "IRAS 12112+0305", "MRK 0231", "WKK 2031", "MRK 0273", "WISEA J143738.49-150019.1", "WISEA J144059.01-370431.8", "IRAS 15250+3609",
          "ARP 220", "2MASX J17232194-0017009", "2MASX J19322229-0400010", "IRAS 19542+1110", "ESO 286-IG 019", "IRAS 22491-1808", "ESO 148-IG 002", "WISEA J233901.29+362108.4"]
print(datatable.keys())

ulirg_ras = []
ulirg_decs = []
ulirg_redshifts = []

for idx, name in enumerate(fixed_names):
    if name in ulirgs:
        ulirg_ras.append(datatable["RA"][idx])
        ulirg_decs.append(datatable["Dec"][idx])
        ulirg_redshifts.append(datatable["Redshift (z)"][idx])
ulirg_ras = np.array(ulirg_ras)
ulirg_decs = np.array(ulirg_decs)
ulirg_redshifts = np.array(ulirg_redshifts)

plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(projection="aitoff"))

gal_l = []
gal_b = []
ulirg_gal_l = []
ulirg_gal_b = []

for idx, ra in enumerate(obj_ras):
    sk = SkyCoord(f"{obj_ras[idx]} {obj_decs[idx]}", unit=(u.deg, u.deg), frame='fk5')
    gal_l.append(sk.galactic.l.wrap_at('180d').radian)
    gal_b.append(sk.galactic.b.radian)

for idx, ra in enumerate(ulirg_ras):
    sk = SkyCoord(f"{ulirg_ras[idx]} {ulirg_decs[idx]}", unit=(u.deg, u.deg), frame='fk5')
    ulirg_gal_l.append(sk.galactic.l.wrap_at('180d').radian)
    ulirg_gal_b.append(sk.galactic.b.radian)

ax.scatter(gal_l, gal_b, marker='o', s=55, color="white", zorder=0)
#cbar = ax.scatter(gal_l, gal_b, c=obj_redshifts, marker='o', s=50, cmap="cool", label="LIRGS", zorder=1)
ax.scatter(gal_l, gal_b, color="cyan", marker='o', s=50, label="LIRGS", zorder=1)
ax.scatter(ulirg_gal_l, ulirg_gal_b, c="gold", s=300, marker='*', label="ULIRGS", zorder=2)#, cmap="cool")
#plt.colorbar(cbar, orientation="horizontal")
ax.set_title("Great Observatories All-sky LIRG Survey (GOALS) Sample", fontsize=24, pad=25)
ax.grid()
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='x', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=18)
plt.legend(loc=(-0.05, -0.14), prop={'size': 16})
plt.savefig("./../diagnostic_plots/goals_aitoff.pdf", dpi=1000)
plt.close()