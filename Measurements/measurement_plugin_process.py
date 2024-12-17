import numpy as np
import matplotlib.pyplot as plt


with np.load('results/shookOnes_wet.npz') as data:

    #measurement parameters
    t = data['t']
    fs = data['fs']
    G_hp = data['G_hp']

    #acquired data
    u_in = data['u_in']
    u_lim = data['u_lim']
    u_hp = data['u_hp']
    x = data['x']

""" PLOT  """
plt.close("all")

fig, ax = plt.subplots(4, sharex=True, layout='constrained') #, sharey=True
ax[0].plot(t, u_in, label='not processed')
ax[0].set(ylabel=r'$U_{in}$ [V]')
ax[0].grid()

ax[1].plot(t, u_lim, label='possibly processed')
ax[1].set(ylabel=r'$U_{lim}$ [V]')
ax[1].grid()

ax[2].plot(t, u_hp, label='Speaker terminal voltage')
ax[2].set(ylabel=r'$U_{hp}$ [V]')
ax[2].grid()

ax[3].plot(t-1e-3, 1000*x)  # compensating displacement sensor delay of 1 ms
#ax[3].plot(t-1e-3, 1000*x_smooth, label='Displacement')
ax[3].set(xlabel='Time [s]', ylabel='Displacement [mm]')
ax[3].grid()

plt.show()