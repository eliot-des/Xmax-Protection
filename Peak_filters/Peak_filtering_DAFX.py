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
G = -25 # Gain (dB)


# model parameter
d = -np.cos(2*np.pi*fc/Fs)
V0 = 10**(G/20)
H0 = V0 -1

c = (np.tan(np.pi*fb/Fs) -V0)/(np.tan(np.pi*fb/Fs) +V0) # cut case


# filters coeff
b0 = 1+H0/2*(1+c)
b1 = d*(1-c)
b2 = -(c+H0/2*(1+c))
a0 = 1
a1 = d*(1-c)
a2 = -c

sos = np.array([[b0, b1, b2, a0, a1, a2]])


# get frequency response
_, H = signal.sosfreqz(sos, worN = freq, fs = Fs)


# Plot it
fig, ax = plt.subplots(2, 1)

ax[0].semilogx(freq, 20*np.log10(np.abs(H)))
ax[1].semilogx(freq, np.rad2deg(np.angle(H)))

ax[0].axhline(y=-3, color='red', linestyle='--', label='-3 dB cutoff')
ax[0].set(xlim = (20, 210), ylim = (-45, 3))
ax[0].set_ylabel('Magnitude [dB ref=1]')
ax[0].set_title(f'Peak filter response. Fc: {fc} Hz - Bandwidth: {fb} Hz')
ax[0].grid(which='both', axis='both')
ax[0].legend(loc='lower right')

ax[1].set(xlabel = 'Frequency [Hz]', ylabel = 'Angle [deg]')
ax[1].set(xlim = (20, 210), ylim = (-90, 90))
ax[1].grid(which='both', axis='both')

plt.show()