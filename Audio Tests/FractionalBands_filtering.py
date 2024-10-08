import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

"""
Filter bandpass bank based on Butterworth response
sos format is used for order higher than 2
"""


Fs = 48000 # sampling freq
freq = np.arange(20, Fs/2, 0.5) # freq tab

fstart = 20 # lower cutoff frequency
fstop = 100 # highest cutoff frequency
df = 10     # band width

order = 6 # order of the filter

band_low = np.arange(fstart, fstop, df) # array containing the low cutoff freq of each bandpass filter
band_high = np.arange(fstart+df, fstop+df, df) # array containing the high cutoff freq of each bandpass filter

# print(band_low)
# print(band_high)
# print(band_low+df/2) # central frequencies


# --------- initialization --------- #
b = np.zeros((2*order+1, len(band_low))) # filter b coeff are stored in a column. They are as many columns as bands
a = np.zeros((2*order+1, len(band_low))) # filter a coeff are stored in a column. They are as many columns as bands
w = np.zeros((len(band_low), len(freq))) # contains the freq at which the freq. resp. will be computed
h = np.zeros((len(band_low), len(freq)), dtype=complex) # line n°i contains the freq. resp. (complex number) of the correspopnding band filter
sos_bank = [np.zeros((order, 6)) for _ in range(len(band_low))]


# --------- compute digital filters coeff --------- #
for i in range(len(band_low)):

    if order<=2:
        b[:,i], a[:,i] = signal.butter(order, [band_low[i], band_high[i]], btype='bandpass', fs=Fs)
        w[i,:], h[i,:] = signal.freqz(b[:,i], a[:,i], worN=freq, fs=Fs)

    else:
        sos = signal.butter(order, [band_low[i], band_high[i]], btype='bandpass', output='sos', fs=Fs)
        sos_bank[i] = sos
        w[i,:], h[i,:] = signal.sosfreqz(sos, worN=freq, fs=Fs)



# --------- plot the freq. resp. --------- #
fig, ax = plt.subplots(2, 1)

for i in range(len(band_low)):

    ax[0].semilogx(w[i,:], 20*np.log10(np.abs(h[i,:])))
    ax[1].semilogx(w[i,:], np.rad2deg(np.angle(h[i,:])))


ax[0].set_ylim([-40, 3])
ax[0].set_xlim([20, 1000])
ax[0].set_ylabel('Magnitude [dB ref=1]')
ax[0].grid(which='both', axis='both')

ax[1].set_xlabel('Frequency (Hz)')
ax[1].set_ylabel('Angle (deg)')
ax[1].set_xlim([20, 1000])
ax[1].grid(which='both', axis='both')

plt.show()