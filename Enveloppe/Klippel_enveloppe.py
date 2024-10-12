import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize
import sys

# Import data and shorten
Fs, v = read('Audio Tests/ShookOnesDisplacement.wav')
t = np.arange(0, len(v)/Fs, 1/Fs)
v = normalize(v)    #normalize the signal to 1

tstart = 1.8
tstop = 2.3

t = t[int(tstart*Fs):int(tstop*Fs)]
v = v[int(tstart*Fs):int(tstop*Fs)]


# Define parameters
Xmax = 0.6
Fr = 72 # resonance frequency used for the displacement filter design


# Compute enveloppe from Klippel trick
A = np.zeros_like(v)

for i in range(len(v)-1):
    A[i] = np.sqrt( v[i]**2 + (1/(2*np.pi*Fr)*(v[i+1]-v[i])*Fs)**2 )


# Design second order low pass FIR filter (with 3 taps) and get frequency response
numtaps = 100  # Second order FIR = 3 taps
cutoff = 100  # Low pass cutoff frequency
a = sig.firwin(numtaps, cutoff, pass_zero='lowpass', fs=Fs)

f = np.geomspace(20, Fs/2, 1000)
_, h = sig.freqz(a, [1], worN=f, fs=Fs)

# Plot frequency response
fig, ax = plt.subplots()
ax.semilogx(f, 20*np.log10(np.abs(h)))
ax.axhline(y=-3, color='red', linestyle='--', label='-3 dB')
ax.set_ylim([-9, 1])
ax.set_title('FIR Low-Pass frequency response')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Amplitude [dB]')
plt.grid(which='both')
plt.show()

# sys.exit()


# Apply the filter to the envelope
A_LP = sig.lfilter(a, [1], A)


# Plot everything
fig, ax = plt.subplots()
ax.plot(t, v, label='Displacement signal')
ax.plot(t, A, label='Enveloppe')
ax.plot(t, A_LP, label='LP filtered enveloppe')
# ax.axhline(y=Xmax, color='red', linestyle='--', label='Xmax treshold')
ax.set_title('Enveloppe follower Klippel')
ax.set_xlabel('Time [sec]')
ax.set_ylabel('Amplitude [a.u]')
ax.legend()
ax.grid()
plt.show()