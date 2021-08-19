import numpy as np
import openmc
import openmc.lib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import cauchy
import csv
import math

plt.rcParams.update({'font.size': 16})

openmc.lib.init()

nuclide = openmc.lib.nuclides['Th232']

# Need to set directory to be same at the file you are accessing for nuclide
nuclide_hdf5 = openmc.data.IncidentNeutron.from_hdf5('/home/jacob/nuclear_data/endf71_unresolved/endf-b7.1-hdf5/neutron/Th232.h5')

# Check which resonance ranges we have stored
# print(vars(nuclide_hdf5._resonances))
# Index manually set for now


unresolved = nuclide_hdf5._resonances._ranges[0]
if not (type(unresolved) == openmc.data.resonance.Unresolved):
    raise TypeError('Resonance object is not "unresolved"')
# print(vars(unresolved)) # <- use to explore data stored in object


indexed_energies = unresolved.energies
energy_idx = 0
energy = indexed_energies[energy_idx]
print(energy)


T = 3000
xs = nuclide.sample_urr_njoy(energy_idx, T, n_sample=1000000,prn_seed=None)
elast_xs = xs.values[:,1]
capture_xs = xs.values[:,3]
print(xs.values)

# Making histogram from cross section samples
hist_elast, bin_edges_elast = np.histogram(elast_xs, bins=10000, density=True)
hist_cap, bin_edges_cap = np.histogram(capture_xs, bins=10000, density=True)

bin_centers_elast = (bin_edges_elast[:-1] + bin_edges_elast[1:])/2
bin_centers_cap = (bin_edges_cap[:-1] +  bin_edges_cap[1:])/2

plt.plot(bin_centers_elast, hist_elast, label= 'OpenMC-Generated PDF')


def cauchyf(x, loc, scale):
    return 	cauchy.pdf(x, loc, scale)


optimalparameters = cauchy.fit(elast_xs, loc= 0, scale= 1.1017170250507625e-06)


# Custom Sampling
sampled_list = []
mean = optimalparameters[0]
scale = optimalparameters[1]

#def sampler():
#	value = np.random.uniform(0, 1)
#	return scale * math.tan(value*(math.pi/2 - math.atan(-mean/scale)) + math.atan(-mean/scale)) + mean

cauchy_obj = cauchy.rvs(loc=mean, scale=scale, size=1000000)

#for i in range(0, 100000):
#	sampled_list.append(sampler())

corrected_sampled_list = cauchy_obj#[]
#for el in sampled_list:
#	if el <= max(bin_centers_cap):
#		corrected_sampled_list.append(el)

sampled_hist, bin_edges_sampled = np.histogram(corrected_sampled_list, bins=bin_edges_elast, density = True)
print(scale)
print(mean)
#print(bin_edges_sampled)
#print(bin_edges_cap)
bin_centers_sampled = (bin_edges_sampled[:-1] + bin_edges_sampled[1:])/2
plt.plot(bin_centers_sampled, sampled_hist, label = 'PDF Sampled From Cauchy Fit')
plt.plot(bin_centers_elast, sampled_hist - hist_elast, label='PDF Minus Cauchy Fit')

plt.title('Cauchy Fit of Th232 OpenMC-Generated Capture Cross Section PDF, 4000 eV, 0K')
plt.ylabel('Probability Density')
plt.xlabel('Th232 Cross Section Values')
plt.legend()
plt.show()
