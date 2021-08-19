import numpy as np
import openmc
import openmc.lib
from openmc.data.urr import ProbabilityTables
from matplotlib import pyplot as plt
import os
import pickle

# Imports a CDF template, and number of bins. constructs 
def makebinprobs(n_samples, n_bins):
	# ensure 10 values in 1/200
	if n_bins < 10:
		raise ValueError('n_bins must be 10 or greater, currently: {}'.format(n_bins))

	fraction = np.array([1/200, 1/40, 1/10, 1/4, 1/2])
	inner = np.ones(n_bins - 10)
	full_fraction = np.hstack([fraction, inner, fraction[::-1]])
	#print(full_fraction)
	norm_fraction = full_fraction/np.sum(full_fraction)
	#print(norm_fraction)
	#if n_samples*norm_fraction[0] < 10:
	#	raise ValueError('Edge bins need at minimum 10 values. Current: {}'.format(int(n_samples*norm_fraction[0])))
	samples_left = n_samples
	bin_size_array = np.zeros(n_bins)
	for i in range(n_bins):
		bin_size = np.floor(n_samples*norm_fraction[i])
		samples_left -= bin_size
		bin_size_array[i] = bin_size
		#print('Samples left: ', samples_left)
	#print(sum(bin_size_array))
	middleidx = n_bins//2
	bin_size_array[middleidx] +=samples_left 
	return bin_size_array




# Samples over each energy and temperature. Note:doing total xs
def samplingroutine(n_bins, energy_idx, temp, n_samples):

	bin_size_array = makebinprobs(n_samples, n_bins)


	xs = nuclide.sample_urr_njoy(energy_idx, temp, n_sample=n_samples, prn_seed=None)
	#print(xs)
	xs = xs.values
	xs = xs[xs[:,0].argsort()]
	#print(bin_size_array)	
	#print(xs)
	def make_xs_bins(xs_array, bin_size_array):

		numbins = len(bin_size_array)
		temporary_array = np.zeros(numbins)
		sum_array = np.zeros(numbins)

		bin_size_array = bin_size_array[::-1]

		for i in range(len(bin_size_array)):
			#print(bin_size_array[i])
			temporary_array = xs_array[-int(bin_size_array[i]):]
			#print(len(temporary_array))
		
			xs_array = xs_array[:-int(bin_size_array[i])]
			#print(len(xs_array))
			#print(xs_array)

			#print(temporary_array)

			for element in temporary_array:
				sum_array[i] += element[0]


		for i in range(len(sum_array)):
			sum_array[i] /= bin_size_array[i]

		return sum_array[::-1]

	return make_xs_bins(xs, bin_size_array)


# returns average, high, and low for each bin in imported array
def analyzebins(arr, bin_size_array):
	average_array = np.zeros(bin_size_array)
	for i in range(len(average_array)):
		average_array[i] = arr / bin_size_array[i]

	return average_array

# 40 total temperatures to start


def get_mesh_array(temp_array, energy_array, numbins, n_samples):


	results_array = np.zeros([len(energy_array), len(temp_array), numbins])

	for e, energy in enumerate(energy_array):
		for t, temp in enumerate(temp_array):
			results_array[e][t][:] = samplingroutine(numbins, e, temp, n_samples)

	return results_array

# Definitions
os.environ['OPENMC_CROSS_SECTIONS'] = '/home/jacob/nuclear_data/endf71_unresolved/endf-b7.1-hdf5/cross_sections.xml'

openmc.lib.init()
isotope = 'Pu239'
nuclide = openmc.lib.nuclides[isotope]
nuclide_hdf5 = openmc.data.IncidentNeutron.from_hdf5('/home/jacob/nuclear_data/endf71_unresolved/endf-b7.1-hdf5/neutron/{}.h5'.format(isotope))
unresolved = nuclide_hdf5._resonances._ranges[0]



energy_array = unresolved.energies[:40]

log_temp_array = np.logspace(0, 6, 40)

x, y = np.meshgrid(energy_array, log_temp_array)

n_bins = 40
n_samples = 25000

binarr = makebinprobs(n_samples, n_bins)

array = get_mesh_array(log_temp_array, energy_array, n_bins, n_samples)

data = {'n_bins': n_bins, 
	'n_samples': n_samples,
	'log_temp_array': log_temp_array,
	'energy_array': energy_array,
	'xstable': array,
	'isotope': isotope
}

filename = '{}_xstable.pkl'.format(isotope)
with open(filename, 'wb') as f:
    pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

	