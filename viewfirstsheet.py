import numpy as np
import pickle
from matplotlib import pyplot as plt

filename = 'Pu239_xstable.pkl'
with open(filename, 'rb') as f:
    array_loaded = pickle.load(f)


z = array_loaded['xstable'][:,:,-1]
x, y = np.meshgrid(array_loaded['energy_array'], array_loaded['log_temp_array'])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z)
ax.set_title("First Sheet Pu239 Cross Section Surface")
ax.set_xlabel("Energy Value, eV")
ax.set_ylabel("Log10 Temperature Value, K")
ax.set_zlabel("First Bin Value")
plt.show()