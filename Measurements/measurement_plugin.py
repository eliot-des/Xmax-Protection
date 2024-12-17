import numpy as np
import matplotlib.pyplot as plt
from functions.measurement_NI import measurement_NI

""" Measurement Parameters """

music = 'shookOnes'
threshold = 0.1  # [V]




""" Parameters """
Dev = 'Dev1'  # name of the NI device
fs = 48000     # sampling frequency (very low to avoid huge data files)

""" Sensitivities """
voltage_sensitivity = 1  # [V/V]
current_sensitivity = 1  # [A/V]
#displacement_sensitivity = 2e-3  # [m/V]
#displacement_DC_offset = 2.5  # offset of the displacement sensor [V]


displacement_sensitivity = 1280 # [um/V]
displacement_sensitivity*=1e-6  # [m/V]

""" Prepare a signal of 5 minutes (just zeros) """
duration = 5  # [s]
t = np.arange(0, duration, 1/fs)  # time axis [s]

output = np.zeros(len(t))


""" Measurement  """
y = measurement_NI(output, fs, Dev)

u_in  = np.array(y[0])     #[V]
u_lim = np.array(y[1])     #[V]
u_hp  = np.array(y[2])     #[V]
x     = np.array(y[3])        
x *= displacement_sensitivity  #[m]


print('Max sortie carte son:', np.max(u_in))
print('Max HP:', np.max(u_hp))
G_amp = np.max(u_hp)/np.max(u_in)
G_hp = np.max(u_hp)


""" SAVE  """
np.savez('results/shookOnes_wet.npz', u_in=u_in, u_lim=u_lim, u_hp = u_hp, x=x, x_smooth=x_smooth, t=t, fs=fs, G_hp = G_hp)

""" PLOT  """
plt.close("all")

fig, ax = plt.subplots(4, sharex=True, sharey=True, layout='constrained') #, sharey=True
ax[0].plot(t, u_in, label='not processed')
ax[0].set(ylabel=r'$U_{in}$ [V]')
ax[0].grid()

ax[1].plot(t, u_lim, label='possibly processed')
ax[1].set(ylabel=r'$U_{lim}$ [V]')
ax[1].grid()

ax[2].plot(t, u_hp, label='Speaker terminal voltage')
ax[2].plot(t, u_lim*G_amp,alpha=0.5, label='Speaker terminal voltage')
ax[2].set(ylabel=r'$U_{hp}$ [V]')
ax[2].grid()

ax[3].plot(t-1e-3, 1000*x)  # compensating displacement sensor delay of 1 ms
#ax[3].plot(t-1e-3, 1000*x_smooth, label='Displacement')
ax[3].set(xlabel='Time [s]', ylabel='Displacement [mm]')
ax[3].grid()

fig, ax = plt.subplots()
ax.plot(t-1e-3, 1000*x)  # compensating displacement sensor delay of 1 ms
ax.set(xlabel='Time [s]', ylabel='Displacement [mm]')
ax.grid()

plt.show()