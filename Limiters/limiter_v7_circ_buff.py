import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write

import os
import sys
# insert root directory into python module search path
fpath = os.path.join(os.getcwd(), 'Modules')
sys.path.append(fpath)

from audio import normalize

#set tex font for plots
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 18,
    "legend.fontsize": 16})

class DelayLine:

    def __init__(self, maxSamplesNbr, initValues = 0):
        self.maxSamplesNbr = maxSamplesNbr
        self.buffer = np.ones(maxSamplesNbr) * initValues
        self.writeIndex = 0

    def write(self, x):
        self.buffer[self.writeIndex] = x
        self.writeIndex = (self.writeIndex + 1) % self.maxSamplesNbr
    
    def read(self, sampleDelay):
        readIndex = (self.writeIndex - sampleDelay) % self.maxSamplesNbr
        return self.buffer[readIndex]

    def getBuffer(self):
        return self.buffer


t_max = 0.1
t_step_start = t_max/3
t_step_stop  = 2*t_max/3
step_amp = 0.2

Fs = 44100
t = np.arange(0, t_max, 1/Fs)
n = np.arange(0, len(t))

f = 20/t_step_start
x = np.sin(2*np.pi*f*t)
x = x*(step_amp + (1-step_amp)*(n>=t_step_start*Fs)*(n<=t_step_stop*Fs))

#Test with a dirac

t_dirac = 0.01

x[0:int(t_step_start*Fs/2)] = 0
x[int(t_dirac*Fs)] = 1

'''
#=============================================================

music = 'Sacrifice1'
Fs, x = read(f'Audio/{music}.wav') 
x = x[:,0]          #select only one channel
x = normalize(x)    #normalize the signal to 1
G = 1               #gain of the amplifier -> Max tension in volts
x*=G             #tension in volts

t = np.arange(0, len(x)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 3            #duration time in seconds

x = x[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]
'''
#=============================================================
#=============================================================

# Define Limiter Threshold
thres = 0.4

# Envelope estimation parameters
attack_time   = 0.002
hold_time     = 0#0.004
release_time  = 0.003#0.003

release_coeff = 1 - np.exp(-2.2/(Fs*release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
x_lim = np.zeros_like(x)
x_g   = np.ones_like(x)
c     = np.ones_like(x)
m     = np.ones_like(x)
g     = np.ones_like(x)



gc = 0
C  = 1
G = 1

#circular buffers:
delayLine = DelayLine(N_attack)
gainLine = DelayLine(N_attack + N_hold, 1)
averageLine = DelayLine(N_attack, 1)

#=============================================================
# Implement the limiter in a for loop

for n in range(1, len(x)):

    delayLine.write(x[n])

    gC = np.minimum(1, thres/np.abs(x[n]))
    x_g[n] = gC
    gainLine.write(gC)

    minGain = np.min(gainLine.getBuffer())
    m[n] = minGain

    
    C = np.minimum(minGain, (1-release_coeff)*C + release_coeff*minGain)
    c[n] = C

    CPrev = averageLine.read(N_attack)
    averageLine.write(C)

    G = G + (1/N_attack) * (C - CPrev) #don't seems to work if we have time varying control parameters...
    g[n] = G
    
    x_lim[n] =  G * delayLine.read(N_attack)

#=============================================================

fig, ax = plt.subplots(3, 1, figsize=(12,8), sharex=True, layout='constrained')

#fig.suptitle(f'Sinusoidal signal of frequency {f:.2f} Hz with amplitude modulation.')

ax[0].plot(t, x, label=r'$x[n]$')
ax[0].plot(t, np.roll(x_lim, -N_attack),'k', label=r'$x_{lim}[n+N_{attack}]$')
ax[0].plot(t, thres*np.ones_like(t), 'r--')
ax[0].plot(t, -thres*np.ones_like(t), 'r--')
ax[0].legend(loc='lower left')

ax[1].plot(t, x_g, label='Gain computer')
ax[1].plot(t, m, '--',label=r'+ Min. filter')
ax[1].plot(t, c, ':', label=r'+ Release')
ax[1].plot(t, g, color='crimson', label=r'+ Average smooth.')
ax[1].legend(loc='center right')

ax[2].plot(t, x_lim, 'k', label=r'$x_{lim}[n]$')
ax[2].plot(t, thres*np.ones_like(t), 'r--')
ax[2].plot(t, -thres*np.ones_like(t), 'r--')
ax[2].set(xlabel='Time [s]')
ax[2].legend(loc='lower right')

ax2twin = ax[2].twinx()
ax2twin.plot(t, g, color='crimson',alpha=0.5, label='Gain function')
ax2twin.set(ylabel='Gain')
ax2twin.legend(loc='upper right')

for i in range(3):
    ax[i].set(ylabel='Amplitude [mm]')
    ax[i].grid()
    ax[i].set_xlim([0, t[-1]])


plt.show()

fig.savefig("Limiter_noHold.pdf", format="pdf", transparent=True)

#=============================================================
#Export the output signal without normalizing it

'''
x_lim = x_lim * 32767
x_lim = x_lim.astype(np.int16)

#write the output signal
write(f'Audio/Limiter/{music}.wav', Fs, x_lim)
'''