import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

# ================================================================
with np.load('Measurements/results/shookOnes_dry.npz') as data:

    #measurement parameters
    t = data['t']
    Fs = int(data['fs'])
    G_hp = data['G_hp']

    #acquired data
    u_in = data['u_in']
    u_lim = data['u_lim']
    u_hp = data['u_hp']
    x = data['x']

# ================================================================
# Load the loudspeaker parameters
with open(f'Dataset_T&S/Peerless_HDSP830860_Klippel.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)

# Low-frequency approximation of X/U(s)
b_xu = np.array([0, 0, Bl/Rec])         
a_xu = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])
bd_xu, ad_xu = sig.bilinear(b_xu, a_xu, fs=Fs)

x_prediction = sig.lfilter(bd_xu, ad_xu, u_hp)
# ================================================================
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

ax[3].plot(t-0.3e-3, 1e3*x, label='Measurements')   #delayed by 0.3 ms -> lazer delay
ax[3].plot(t, 1e3*x_prediction, label=r'Predicted from $U_{hp}$')
ax[3].set(xlabel='Time [s]', ylabel='Displacement [mm]')
ax[3].grid()
ax[3].legend()
plt.show()

'''
fig, ax = plt.subplots()

r_xxpred = np.correlate(x, x_prediction, mode='full')
l = np.arange(-len(r_xxpred)//2, len(r_xxpred)//2)
lmax = l[np.argmax(r_xxpred)]
print(f"Predicted delay: {lmax/Fs} s")
ax.plot(l, r_xxpred)
plt.show()
'''