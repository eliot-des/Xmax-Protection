import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sys

"""
Frequency response of Peak filters 
Filter model based on DAFX book, P64 - P83 PDF
"""

# sampling frequency
Fs = 48000


# frequency tab
# freq = np.geomspace(1, Fs/2, 1000)
freq = np.geomspace(20, 200, 1000)


# filter parameters
fc = 50 # central frequency
fb = 10 # bandwidth
G_list = np.round(np.geomspace(-6, -40, 5)).astype(int)
print(G_list)


# model parameter
d = -np.cos(2*np.pi*fc/Fs)
V0 = 10**(G_list/20)
H0 = V0 -1

c = np.array([(np.tan(np.pi*fb/Fs) -i)/(np.tan(np.pi*fb/Fs) +i) for i in V0]) # cut case


# filters coeff
b0 = np.array([1+j/2*(1+k) for j, k in zip(H0, c)])
b1 = np.array([d*(1-k) for k in c])
b2 = np.array([-(k+j/2*(1+k)) for j, k in zip(H0, c)])
a0 = np.ones(len(G_list))
a1 = b1
a2 = np.array([-k for k in c])

# Concatenate b0, b1, b2 into a matrix B where each row is [b0, b1, b2]
B = np.column_stack([b0, b1, b2])

# Concatenate a0, a1, a2 into a matrix A where each row is [a0, a1, a2]
A = np.column_stack([a0, a1, a2])

# creat non-normalized SOS matrix
sos = np.concatenate([B, A], axis=1)  # already normalized as a0 = 1 for all


H = np.zeros((sos.shape[0], len(freq)), dtype = complex)

for i in range(sos.shape[0]):
    _, H[i] = signal.sosfreqz(sos[i,:], worN = freq, fs = Fs)


#create color map for the bands
colors = plt.cm.viridis(np.linspace(0.1, 0.9, sos.shape[0]))


# Plot it
fig, ax = plt.subplots(2, 1)

for i in range(sos.shape[0]):
    ax[0].semilogx(freq, 20*np.log10(np.abs(H[i])), color = colors[i], label=f'Gain: {G_list[i]} dB')
    ax[1].semilogx(freq, np.rad2deg(np.angle(H[i])), color = colors[i])

ax[0].axhline(y=-3, color='red', linestyle='--', label='-3 dB cutoff')
ax[0].set(xlim = (20, 210), ylim = (-45, 3))
ax[0].set_ylabel('Magnitude [dB ref=1]')
ax[0].set_title(f'DAFX P64 peak filter response. Fc: {fc} Hz - Bandwidth: {fb} Hz')
ax[0].grid(which='both', axis='both')
ax[0].legend(loc='lower right')

ax[1].set(xlabel = 'Frequency [Hz]', ylabel = 'Angle [deg]')
ax[1].set(xlim = (20, 210), ylim = (-90, 90))
ax[1].grid(which='both', axis='both')

# plt.savefig('Figures/Peak_filtering_DAFX.pdf')
plt.show()