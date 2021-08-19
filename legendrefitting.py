import numpy as np
import pickle
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

filename = 'Pu239_xstable.pkl'
with open(filename, 'rb') as f:
    array_loaded = pickle.load(f)

X, Y = np.meshgrid(array_loaded['energy_array'], array_loaded['log_temp_array'])
cross_sections = array_loaded['xstable'][:,:,-1]

# plt.show()

def leg2d(loc, *args):
    # Function callable by optimizer
    vec = np.array(args).flatten()
    shape = int(np.sqrt(len(vec)))
    c = vec.reshape([shape, shape])
    x, y = loc
    # print('input', x, y, c)
    return np.polynomial.legendre.legval2d(x, y, c)
    # return np.polynomial.polynomial.polyval2d(x, y, args)

# Test legendre
# nx = 10
# ny = 10
# x = np.linspace(-1, 1, nx)
# y = np.linspace(-1, 1, ny)
# testX, testY = np.meshgrid(x, y)
# testpoints = np.vstack((testX.ravel(), testY.ravel()))
# # testpoints = [x, y]
# c = np.zeros([3,3])
# c[0,0] = 1
# c[1,1] = 1
# testZ = leg2d(testpoints, c).reshape(ny, nx)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(testX, testY, testZ, cmap='viridis')
# plt.show()

def map_domain(x):
    return 2*(x-x[0])/(x[-1] - x[0])-1

# Map domain to [-1, 1] for legendre fit
legx = map_domain(array_loaded['energy_array'])
legy = map_domain(array_loaded['log_temp_array'])
legX, legY = np.meshgrid(legx, legy)

# Initialize guess of parameters
order = 30
c_guess = np.ones(order*order)

# Fit data using scipy.optimize
xdata = np.vstack((legX.ravel(), legY.ravel()))
c_opt, c_cov = curve_fit(leg2d, xdata, cross_sections.ravel(), 
                         c_guess, maxfev = 100000)

print('Optimal c matrix\n', c_opt)
c_rel_error = list(c_cov[i,i]/c_opt[i] for i in range(order**2))
c_rel_error = np.array(c_rel_error).reshape(order,order)
with np.printoptions(precision=1):
    print('Relative Error on c\n', c_rel_error)

# Plotting
# Custom legend
custom_lines = [Line2D([0], [0], color='b', lw=4),
                Line2D([0], [0], color='r', lw=4)]

fig = plt.figure(figsize=(9, 4))
ax = fig.add_subplot(121, projection='3d')
fit_xs = leg2d(xdata, c_opt).reshape([len(array_loaded['log_temp_array']), len(array_loaded['energy_array'])])

plot_x = X
plot_y = Y
# ax.plot_surface(X, Y, cross_sections, cmap='Blues')
# ax.plot_surface(X, Y, leg_xs, cmap='Reds')
ax.plot_surface(plot_x, plot_y, cross_sections, cmap='Blues')
print('legX', legX.shape, 'legY', legY.shape, 'fit_xs', fit_xs.shape)
ax.plot_surface(plot_x, plot_y, fit_xs, cmap='Reds')
ax.legend(custom_lines, ['Example Data', 'Fit'])
ax.set_xlabel('E')
ax.set_ylabel('log(T)')
ax.set_zlabel('Cross Section (barns)')

ax = fig.add_subplot(122)
contour = ax.contourf(plot_x, plot_y, (fit_xs - cross_sections)/cross_sections*100, 
                    cmap = "seismic")
cbar = fig.colorbar(contour)
cbar.set_label('% Relative Difference')
ax.set_title('Relative Difference')
ax.set_xlabel('E')
ax.set_ylabel('log(T)')

fig.tight_layout()
plt.show()