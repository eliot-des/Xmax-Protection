import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write

import os
import sys
# insert root directory into python module search path
fpath = os.path.join(os.getcwd(), 'Modules')
sys.path.append(fpath)

from audio import normalize

print_test()
print_test2()

#=============================================================

#Implement a basic hard-clipping, for which the side chain is:

#-take the absolute value of the input signal -> abs_x
#-apply the gain computer function to the absolute value -> x_g
#-apply a minimum filter to the gain computer output -> c
#   -this minimum filter is of lenght N_attack + N_hold !


#-apply an exponential release to the minimum filter output -> c
#-apply FIR filter to this last output with a delay-> g

#=============================================================
#create a sinus signal modulated in amplitude by a step signal

t_max = 0.1
t_step_start = t_max/3
t_step_stop  = 2*t_max/3
step_amp = 0.2

Fs = 44100
t = np.arange(0, t_max, 1/Fs)
n = np.arange(0, len(t))

f = 2/t_step_start
x = 1#np.sin(2*np.pi*f*t)
x = x*(step_amp + (1-step_amp)*(n>=t_step_start*Fs)*(n<=t_step_stop*Fs))

#Test with a dirac
x = np.zeros_like(t)
x[int(t_step_start*Fs)] = 1


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

#=============================================================
#=============================================================

# Define Limiter Threshold
thres = 0.5

# Envelope estimation parameters
attack_time   = 0.002
hold_time     = 0.002
release_time  = 0.0085

release_coeff = 1 - np.exp(-2.2/(Fs*release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
x_lim = np.zeros_like(x)
x_g   = np.ones_like(x)
c     = np.ones_like(x)
g     = np.ones_like(x)


#=============================================================
# Implement the limiter in a for loop

#with a rectangular FIR filter
fir_coeffs = np.ones(N_attack) / N_attack

#with a hanning FIR filter
#fir_coeffs = np.hanning(N_attack)
#fir_coeffs = fir_coeffs / np.sum(fir_coeffs)

for n in range(N_attack+N_hold, len(x)):
    # Calculate the absolute value of the input signal
    abs_x = np.abs(x[n])

    # Apply the gain computer function to the absolute value
    x_g[n] = np.minimum(1, thres/abs_x)
    #x_g[n] =  gain_computer_compressor(abs_x, thres, 1.5)

    # Apply the minimum filter to the gain computer output
    c[n] = np.min(x_g[n-(N_attack+N_hold):n+1])

    # Apply exponential release to the minimum filter output
    c[n] = np.minimum(c[n], (1-release_coeff)*c[n-1] + release_coeff*c[n]) #recursive implementation

    # Apply the averaging filter to the minimum filter output
    if n >= N_attack:
        #g[n] = np.dot(fir_coeffs, c[n-N_attack+1:n+1])
        g[n] = g[n-1] + 1/N_attack * (c[n] - c[n-N_attack]) #recursive implementation
    else:
        g[n] = c[n]  # No filtering for the initial samples

    # Apply the gain function to the delayed gain computer output
    x_lim[n] = x[n-N_attack]* g[n]

#=============================================================

fig, ax = plt.subplots(3, 1, sharex=True)

#fig.suptitle(f'Sinusoidal signal of frequency {f:.2f} Hz with amplitude modulation.')

ax[0].plot(t, x, label=r'$x[n]$')
ax[0].plot(t, np.roll(x_lim, -N_attack),'k', label=r'$x_{lim}[n]$')
ax[0].plot(t, thres*np.ones_like(t), 'r--', label='Threshold')


ax[1].plot(t, x_g, label='Gain computer output')
ax[1].plot(t, c, label=r'$c[n]$')
ax[1].plot(t, g,color='crimson', label=r'$g[n]$')

ax[2].plot(t, x_lim, label=r'$x_{lim}[n]$')
ax[2].plot(t, thres*np.ones_like(t), 'r--', label='Threshold')
ax[2].set(xlabel='Time [s]')

ax2twin = ax[2].twinx()
ax2twin.plot(t, g, 'g',alpha=0.2, label='Gain function')
ax2twin.set(ylabel='Gain')
ax2twin.legend(loc='upper right')

for i in range(3):
    ax[i].set(ylabel='Amplitude')
    ax[i].legend()
    ax[i].grid()
    ax[i].set_xlim([0, t[-1]])



plt.show()



#=============================================================
#Export the output signal without normalizing it

'''
x_lim = x_lim * 32767
x_lim = x_lim.astype(np.int16)

#write the output signal
write(f'Audio/Limiter/{music}.wav', Fs, x_lim)
'''