import os 
import scipy.interpolate as interpolate 
import numpy as np
import matplotlib.pyplot as plt 

# plt.style.use('bmh')


# Debugging. 
os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf')


dft_energies = np.loadtxt('dft_energies.csv', skiprows=1, delimiter=',')
gw_energies = np.loadtxt('gw_bse_energies.csv', skiprows=1, delimiter=',')

# dft_interp = interpolate.interp1d(dft_energies[:, 0], dft_energies[:, 1], kind='quadratic')
dft_interp = np.poly1d(np.polyfit(dft_energies[:, 0], dft_energies[:, 1], deg=2))
# gw_interp = interpolate.interp1d(gw_energies[:, 0], gw_energies[:, 1], kind='quadratic')
gw_interp = np.poly1d(np.polyfit(gw_energies[:, 0], gw_energies[:, 1], deg=2))

position = 1.15
delta = 0.01
dft_force = - (dft_interp(position + delta) - dft_interp(position))/delta
gw_force = - (gw_interp(position + delta) - gw_interp(position))/delta
print(f'DFT energy at {position} A is: {dft_interp(position)} eV')
print(f'GW-BSE energy at {position} A is: {gw_interp(position) - dft_interp(position)} eV')
print(f'Total energy at {position} A is: {gw_interp(position)} eV')
print(f'DFT force at {position} A is: {dft_force} eV/A')
print(f'GW-BSE force at {position} A is: {gw_force - dft_force} eV/A')
print(f'Total force at {position} A is: {gw_force} eV/A')

# plot gw-bse energies and forces. 
dist_range = np.linspace(1.1, 1.35, 50)
fig = plt.figure(figsize=(12, 6))

ax = fig.add_subplot(1, 2, 1)
ax.plot(dist_range, gw_interp(dist_range), label='Total')
ax.plot(dist_range, dft_interp(dist_range), label='DFT')
ax.plot(dist_range, gw_interp(dist_range) - dft_interp(dist_range), label='GW-BSE alone')
ax.set_title('Energies')
ax.set_xlabel('Bond length (A)')
ax.set_ylabel('Energy (eV)')
ax.legend()
ax.grid()

ax = fig.add_subplot(1, 2, 2)
dft_forces = - (dft_interp(dist_range + delta) - dft_interp(dist_range))/delta 
total_forces = - (gw_interp(dist_range + delta) - gw_interp(dist_range))/delta 
ax.plot(dist_range, total_forces, label='total')
ax.plot(dist_range, dft_forces, label='DFT')
ax.plot(dist_range, total_forces - dft_forces, label='GW-BSE forces alone')
ax.set_title('Forces')
ax.set_xlabel('Bond length (A)')
ax.set_ylabel('Force (eV/A)')
ax.legend()
ax.grid()


plt.show()