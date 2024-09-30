import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sys

"""
Frequency response of Peak filters 
Filter model based on DAFX book, P64 - P83 PDF
References [Whi86, RM87, Dut89a, HB93, Bri94, -Orf96-, -Orf97-, Zol05]
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
w0 = 2*np.pi*fc/Fs
Q = 1 #fc/fb
alpha = np.sin(w0)/(2*Q)
A_dBlist = 10**(G_list/40)
G_fix = -3
A_dB_cst = 10**(G_fix/40)

# filters coeff
b0 = np.array([1+alpha*i for i in A_dBlist])
b1 = -2*np.cos(w0)*np.ones(len(A_dBlist))
b2 = np.array([1-alpha*i for i in A_dBlist])
a0 = (1+alpha/A_dB_cst)*np.ones(len(A_dBlist))
a1 = b1
a2 = (1-alpha/A_dB_cst)**np.ones(len(A_dBlist))

# Concatenate b0, b1, b2 into a matrix B where each row is [b0, b1, b2]
B = np.column_stack([b0, b1, b2])

# Concatenate a0, a1, a2 into a matrix A where each row is [a0, a1, a2]
A = np.column_stack([a0, a1, a2])

# creat non-normalized SOS matrix
sos = np.concatenate([B, A], axis=1)  # already normalized as a0 = 1 for all

for i in range(sos.shape[0]):   # normalization of each raw by its a0 coefficient
    sos[i,:] /= sos[i,3]        # sos[i,3] is the a0 coeff of each filter


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
ax[0].set_title(f'Constant denum gain peak filter response.\n Fc: {fc} Hz - Q: {np.round(Q,2)} - Cst Gain: {np.round(G_fix,2)} dB')
ax[0].grid(which='both', axis='both')
ax[0].legend(loc='lower right')

ax[1].set(xlabel = 'Frequency [Hz]', ylabel = 'Angle [deg]')
ax[1].set(xlim = (20, 210), ylim = (-90, 90))
ax[1].grid(which='both', axis='both')

# plt.savefig('Figures/Peak_filtering_Novak_bis.pdf')
plt.show()