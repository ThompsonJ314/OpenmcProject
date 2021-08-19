import numpy as np
import openmc
import openmc.lib
from openmc.data.urr import ProbabilityTables
from matplotlib import pyplot as plt
import os
os.environ['OPENMC_CROSS_SECTIONS'] = '/home/jacob/nuclear_data/endf71_unresolved/endf-b7.1-hdf5/cross_sections.xml'

openmc.lib.init()
nuclide = openmc.lib.nuclides['Pu239']

nuclide_hdf5 = openmc.data.IncidentNeutron.from_hdf5('/home/jacob/nuclear_data/endf71_unresolved/endf-b7.1-hdf5/neutron/Pu239.h5')

urr = nuclide_hdf5.urr['294K']

unresolved = nuclide_hdf5._resonances._ranges[0]

if not (type(unresolved) == openmc.data.resonance.Unresolved):
    raise TypeError('Resonance object is not "unresolved"')



# Main helper function 
def make_averages_array_from_xs(temp, energy_idx, numbins, numsamples):
	xs = nuclide.sample_urr_njoy(energy_idx, temp, n_sample=numsamples, prn_seed=None)
	xs = xs.values
	xs = xs[xs[:,0].argsort()]
	def make_xs_bins(xs_array):

		temporary_array = np.zeros(numbins)
		sum_array = np.zeros(numbins)
		return_array = np.zeros(numbins)
		samplesperbin = int(numsamples/numbins)

		for i in range(numbins):
			temporary_array = xs_array[-samplesperbin:]
		
			xs_array = xs_array[:-samplesperbin]

			for element in temporary_array:
				sum_array[i]+=element

		for i in range(sum_array.size):
			return_array[i] = sum_array[i] / (samplesperbin)

		return return_array[::-1]

	cumulativeprob = np.zeros(numbins)

	for i in range(cumulativeprob.size):
		cumulativeprob[i] = (i+1)*(1/numbins)

	table = np.zeros((6, numbins))
	table[0] = cumulativeprob

	for i in range(1, 5):
		table[i] = make_xs_bins(xs[:,i-1])

	return table

numbins = 20
numsamples = 30000

energies = unresolved.energies
probtable = np.zeros((energies.size, 6, numbins))


# Heating num
for i in range(energies.size):
	probtable[i] = make_averages_array_from_xs(294, i, numbins, numsamples)
	probtable[i][5] = urr.table[i][5]


probtable = np.array(probtable)
#print(probtable)

# Accesing 1 array in first 2D slice of 3D array
totalarray = probtable[0,1,:]
elasticarray = probtable[0,2,:]
fissionarray = probtable[0,3,:]
capturearray = probtable[0,4,:]
cdf = probtable[0,0,:]

totalarray2 = urr.table[0,1,:]
elasticarray2 = urr.table[0,2,:]
fissionarray2 = urr.table[0,3,:]
capturearray2 = urr.table[0,4,:]
cdf2 = urr.table[0,0,:]
#print(capturearray)


# Making Hybrid Binning
energyidx = 0
CDF_at_e_idx = urr.table[0][0]

sampledxs = nuclide.sample_urr_njoy(energyidx, 294, n_sample=numsamples, prn_seed=None)
print(sampledxs)
sampledxs = sampledxs.values
print(sampledxs[:,3])
sampledxs = sampledxs[sampledxs[:,0].argsort()]
totalxs = sampledxs[:,0]
elasticxs = sampledxs[:,1]
fissionxs = sampledxs[:,2]
capturexs = sampledxs[:,3]

prebinarray1 = [0 for i in range(numbins)]
prebinarray2 = [0 for i in range(numbins)]
prebinarray3 = [0 for i in range(numbins)]
prebinarray4 = [0 for i in range(numbins)]


def make_hybrid_array(prebinarray, myxs):

	for i in range(numbins):
		if i == 0:
			probdensity = CDF_at_e_idx[0]
		else:
			probdensity = CDF_at_e_idx[i] - CDF_at_e_idx[i - 1]

		samplenum = int(np.floor(probdensity*numsamples))
		print(samplenum)

		if i != numbins - 1:

			prebinarray[i] = list(myxs[0: samplenum])
			myxs = myxs[samplenum:]

		else:
			prebinarray[i] = list(myxs)

	returnarray = [0 for i in range(len(CDF_at_e_idx))]
	index = -1
	for bin_ in prebinarray:
		index = index+1
		mysum = 0
		samplequant = len(bin_)
		for value_ in bin_:
			mysum = mysum + value_
		avgval = mysum/samplequant
		returnarray[index] = avgval 

	return returnarray


totalarray3 = make_hybrid_array(prebinarray1, totalxs)
elasticarray3 = make_hybrid_array(prebinarray2, elasticxs)
fissionarray3 = make_hybrid_array(prebinarray3, fissionxs)
capturearray3 = make_hybrid_array(prebinarray4, capturexs)





#print(prebinarray)
#print(capturearray)


fig = plt.figure()

ax = fig.add_subplot(221)

ax.plot(np.hstack([0, totalarray]), np.hstack([0, cdf]), drawstyle='steps',label='customcdf')
ax.plot(np.hstack([0, totalarray2]), np.hstack([0, cdf2]), drawstyle='steps',label='uneditedcdf')
ax.plot(np.hstack([0, totalarray3]), np.hstack([0, CDF_at_e_idx]), drawstyle='steps',label='hybrid')
ax.legend()
ax.set_xlabel("Average Total XS Value")
ax.set_ylabel("CPD Value")
ax.set_title("Total XS")


ax1 = fig.add_subplot(222)

ax1.plot(np.hstack([0, elasticarray]), np.hstack([0, cdf]), drawstyle='steps',label='customcdf')
ax1.plot(np.hstack([0, elasticarray2]), np.hstack([0, cdf2]), drawstyle='steps',label='uneditedcdf')
ax1.plot(np.hstack([0, elasticarray3]), np.hstack([0, CDF_at_e_idx]), drawstyle='steps',label='hybrid')
ax1.legend()
ax1.set_xlabel("Average Elastic XS Value")
ax1.set_ylabel("CPD Value")
ax1.set_title("Elastic XS")


ax2 = fig.add_subplot(223)

ax2.plot(np.hstack([0, fissionarray]), np.hstack([0, cdf]), drawstyle='steps',label='customcdf')
ax2.plot(np.hstack([0, fissionarray2]), np.hstack([0, cdf2]), drawstyle='steps',label='uneditedcdf')
ax2.plot(np.hstack([0, fissionarray3]), np.hstack([0, CDF_at_e_idx]), drawstyle='steps',label='hybrid')
ax2.legend()
ax2.set_xlabel("Average Fission XS Value")
ax2.set_ylabel("CPD Value")
ax2.set_title("Fission XS")


ax3 = fig.add_subplot(224)

ax3.plot(np.hstack([0, capturearray]), np.hstack([0, cdf]), drawstyle='steps',label='customcdf')
ax3.plot(np.hstack([0, capturearray2]), np.hstack([0, cdf2]), drawstyle='steps',label='uneditedcdf')
ax3.plot(np.hstack([0, capturearray3]), np.hstack([0, CDF_at_e_idx]), drawstyle='steps',label='hybrid')
ax3.legend()
ax3.set_xlabel("Average Capture XS Value")
ax3.set_ylabel("CPD Value")
ax3.set_title("Capture XS")



#ax = fig.add_subplot(132)

#ax.plot(np.hstack([0, capturearray]), np.hstack([0, cdf]), drawstyle='steps',label='customcdf')
#ax.plot(np.hstack([0, capturearray2]), np.hstack([0, cdf2]), drawstyle='steps',label='uneditedcdf')
#ax.legend()

plt.show() 







