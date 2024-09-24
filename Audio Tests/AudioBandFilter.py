import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.io.wavfile import read, write  
from lib import *

#import/read wavefile
Fs, x = read('Audio Tests/ShookOnesLPFiltered.wav')
# print(x.dtype)

f  = np.geomspace(10, Fs/2, 1000)
fc = [30, 60, 90, 120, 200]


# sos_lwrl = multibandLinkwitzFilters(8, fc, Fs)
sos_lwrl = multibandEllipticFilters(8, fc, Fs, 1, 60)

x_filt = np.zeros((len(sos_lwrl), len(x)))

for i in range(len(sos_lwrl)):
    x_filt[i] = sig.sosfilt(sos_lwrl[i], x)

    # write(f'Audio Tests/Linkwitz/ShookOnesLPFiltered_bandNo{i}.wav', Fs, x_filt[i].astype(np.int16))
    write(f'Audio Tests/Elliptic/ShookOnesLPFiltered_bandNo{i}.wav', Fs, x_filt[i].astype(np.int16))


'''
filtered_data = filtered_data/np.max(np.abs(filtered_data))
filtered_data = filtered_data*2**15
write('Audio Tests/ShookOnesFiltered.wav', Fs, filtered_data.astype(np.int16))

H_lwrl = np.zeros((len(sos_lwrl), len(f)), dtype = complex)
for i in range(len(sos_lwrl)):
    _, H_lwrl[i] = sig.sosfreqz(sos_lwrl[i], worN = f, fs = Fs)

fig, ax = plt.subplots(2, 1, sharex = True)
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(sos_lwrl)))

for i in range(len(sos_lwrl)):
    ax[0].semilogx(f, 20*np.log10(np.abs(H_lwrl[i])), label = f'Band {i+1}', color = colors[i])
    ax[1].semilogx(f, np.angle(H_lwrl[i]), label = f'Band {i+1}', color = colors[i])

ax[0].semilogx(f, 20*np.log10(np.abs(H_tot)), label = 'Total', color = 'black')

ax[0].legend()
ax[0].grid(which='both')
ax[1].grid(which='both')
ax[0].set_title('Modulus and phase response of the filters')
ax[0].set(xlabel = 'Frequency [Hz]', ylabel = 'Magnitude [dB]')
ax[0].set(xlim = (10, Fs/2), ylim = (-40, 10))
ax[1].set(xlim = (10, Fs/2), ylim = (-np.pi, np.pi))
plt.show()
'''