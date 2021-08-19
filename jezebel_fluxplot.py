import openmc
import numpy as np
import os
import matplotlib.pyplot as plt

os.environ['OPENMC_CROSS_SECTIONS'] = '/home/icmeyer/nuclear_data/endf71_unresolved/endf-b7.1-hdf5/cross_sections.xml'
print(os.environ['OPENMC_CROSS_SECTIONS'])

sp = openmc.StatePoint(f'statepoint.h5')
keff = sp.k_combined
print(f'Final k-effective = {keff}')

# Plot flux tally
tally = sp.get_tally(name='flux')
df = tally.get_pandas_dataframe()
print(tally)
print(df)
elow = df['energy low [eV]'].values
ehigh = df['energy high [eV]'].values
ebins = np.hstack([elow[0], ehigh])
emids = (elow + ehigh)/2
flux_mean = df['mean'].values
flux_sd = df['std. dev.'].values

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel('Flux (particle cm/source particle)')
ax.set_xlabel('Energy (eV)')
ax.set_xscale('log')
ax.set_yscale('log')
flux_plot = np.hstack([flux_mean[0], flux_mean])
# plot mean values
ax.step(ebins, flux_plot, label='flux', color='b')
# plot error bars
ax.errorbar(emids, flux_mean, yerr=flux_sd, fmt='none', capsize=2, color='b')
ax.set_title('Flux for Jezebel with 200 batches, 10000 particles')


sp = openmc.StatePoint(f'statepointmodified.h5')
keff = sp.k_combined
print(f'Final k-effective = {keff}')

# Plot flux tally
tally = sp.get_tally(name='flux')
df = tally.get_pandas_dataframe()
print(tally)
print(df)
elow = df['energy low [eV]'].values
ehigh = df['energy high [eV]'].values
ebins = np.hstack([elow[0], ehigh])
emids = (elow + ehigh)/2
flux_mean = df['mean'].values
flux_sd = df['std. dev.'].values

ax.step(ebins, flux_plot, label='fluxmodified', color='r')
ax.errorbar(emids, flux_mean, yerr=flux_sd, fmt='none', capsize=2, color='r')

plt.show()


