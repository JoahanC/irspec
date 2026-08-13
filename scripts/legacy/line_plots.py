import numpy as np 
import matplotlib.pyplot as plt 

labels_b = ["[MgV]", "[NeV]", "[SIV]", "[NeIII]", "[NeII]"]
vel_disp_b = [555.4449650479809, 541.4839559958773, 530.7871495637856, 265.3199787968706, 211.70727036173867]
rel_vel_b = [215.9705904, -165.66828499, -223.26308725, -116.9985789, -18.95026048]
ip = [141.27, 126.21, 47.22, 63.45, 40.9]
labels_n = ["[ArII]", "[ArIII]", "[NeII]", "[NeIII]", "[SIII]", "[SIV]", "[NeV]", "[NeV]_14", "[NeVI]", "[OIV]", "[MgIV]", "[MgV]"]
ip_n = [27.63, 40.74, 40.9, 63.45, 34.79, 47.22, 126.21, 126.21, 157.93, 77.41, 109.24, 141.27]
rel_vel_n = [1.30210777, -27.51717465, 2.03578897, -25.40669897, -21.95167275, -90.90091679, -52.98816545, -77.47453927, -28.88886464, np.nan, -52.3081531, -59.21960783]
vel_disp_n = [101.5792134971385, 118.69224129967027, 84.99081067720637, 90.87363345950646, 131.0675089418765, 131.0529027887501, 159.27827628022857, 133.78381823823855, 160.81034967266413, np.nan, 142.841525598557, 163.36813666475675]

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)

for idx, line in enumerate(labels_b):
    ax.scatter(ip[idx], vel_disp_b[idx], label=line, s=100)
ax.set_xlabel("Ionization Potential [eV]")
ax.set_ylabel("Velocity Dispersion [km/s]")
ax.set_title("Broad Components")
ax.legend()
plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)

for idx, line in enumerate(labels_n):
    ax.scatter(ip_n[idx], vel_disp_n[idx], label=line, s=100)
ax.set_xlabel("Ionization Potential [eV]")
ax.set_ylabel("Velocity Dispersion [km/s]")
ax.set_title("Narrow Components")
ax.legend()
plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)
for idx, line in enumerate(labels_b):
    ax.scatter(ip[idx], abs(rel_vel_b[idx]), label=line, s=100)
ax.set_xlabel("Ionization Potential [eV]")
ax.set_ylabel("Relative Velocity [km/s]")
ax.set_title("Broad Components")
ax.legend()
plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)
for idx, line in enumerate(labels_n):
    ax.scatter(ip_n[idx], abs(rel_vel_n[idx]), label=line, s=100)
ax.set_xlabel("Ionization Potential [eV]")
ax.set_ylabel("Relative Velocity [km/s]")
ax.set_title("Narrow Components")
ax.legend()
plt.show()