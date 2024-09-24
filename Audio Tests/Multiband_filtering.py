import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.io.wavfile import read, write  
from lib import *
Fs = 48000
#import/read wavefile

Fs, data = read('Audio Tests/ShookOnes.wav')


A = 1                                #gain of the amplifier
u = A*normalize(data[:, 1])          #tension in volts
t = np.arange(0, len(data)/Fs, 1/Fs)


# Thiele-small parameters
# Data from Peerless HDS-P830860 Datasheet
fs = 72
Rec = 6.4
Lec = 0.278e-3
Qms = 2.08
Qes = 0.725
Qts = 0.54
Mms = 8.88e-3
Cms = 560e-6
Rms = 1/(2*np.pi*fs*Cms*Qms)
Bl  = 5.74

#create filter
#analog coefficients
b = np.array([0, 0, 0, Bl])
a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])

#convert to digital filter
bd, ad = sig.bilinear(b, a, Fs)

#apply filter to signal
x = sig.lfilter(bd, ad, u)    #displacement in meters

fig, ax = plt.subplots(2, 1)

ax[0].plot(t, u)
ax[1].plot(t, x*1e3)
ax[0].set_title('Tension')
ax[1].set_title('Displacement')
ax[0].set_ylabel('Amplitude [V]')
ax[1].set_ylabel('Displacement [mm]')
ax[1].set_xlabel('Time [s]')
ax[0].grid()
ax[1].grid()


f = np.geomspace(1, Fs/2, 1000)

#multi-band processing
#create 5 bands -> LPF, BPF, BPF, BPF, HPF
#fc = [25, 200, 700, 4000]
fc = [20, 50, 100, 200]

sos_lwrl = multibandLinkwitzFilters(4, fc, Fs) 
sos_elli = multibandEllipticFilters(5, fc, Fs, 1, 60)

H_lwrl = np.zeros((len(sos_lwrl), len(f)), dtype = complex)
H_elli = np.zeros((len(sos_elli), len(f)), dtype = complex)

for i in range(len(sos_lwrl)):
    _, H_lwrl[i] = sig.sosfreqz(sos_lwrl[i], worN = f, fs = Fs)
    _, H_elli[i] = sig.sosfreqz(sos_elli[i], worN = f, fs = Fs)



print(sos_elli[0].shape)

H_tot = H_lwrl[0] - H_lwrl[1] + H_lwrl[2] - H_lwrl[3] + H_lwrl[4]

fig, ax = plt.subplots(2, 1, sharex = True)

#create color map for the bands
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(sos_lwrl)))

for i in range(len(sos_lwrl)):
    ax[0].semilogx(f, 20*np.log10(np.abs(H_lwrl[i])), label = f'Band {i+1}', color = colors[i])
    #ax[1].semilogx(f, np.angle(H_lwrl[i]), label = f'Band {i+1}')
    ax[1].semilogx(f, 20*np.log10(np.abs(H_elli[i])),'--', label = f'Band {i+1}', color = colors[i])

#ax[0].semilogx(f, 20*np.log10(np.abs(H_tot)), label = 'Total')
#ax[1].semilogx(f, np.angle(H_tot), label = 'Total')

ax[0].legend()
ax[0].grid(which='both')
ax[1].grid(which='both')
ax[0].set_title('Modulus and phase response of the filters')
ax[0].set(xlabel = 'Frequency [Hz]', ylabel = 'Magnitude [dB]')
ax[0].set(xlim = (10, Fs/2), ylim = (-40, 10))
ax[1].set(xlim = (10, Fs/2), ylim = (-40, 10))
plt.show()



'''
#normalize filtered data
filtered_data = x/np.max(np.abs(x))
filtered_data = filtered_data*2**15
write('Audio Tests/ShookOnesFiltered.wav', Fs, filtered_data.astype(np.int16))
'''