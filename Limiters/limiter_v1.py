import numpy as np
import matplotlib.pyplot as plt
import os
import sys
# insert root directory into python module search path
fpath = os.path.join(os.getcwd(), 'Modules')
sys.path.append(fpath)

from dynamic_trackers import gain_factor_smoothing

#=============================================================

#Implement a basic hard-clipping limiter, without any 
#lookahead delay

#=============================================================
#create a sinus signal modulated in amplitude by a step signal

t_max = 0.1
t_step_start = t_max/3
t_step_stop  = 2*t_max/3
step_amp = 0.2

Fs = 44100
t = np.arange(0, t_max, 1/Fs)
n = np.arange(0, len(t))

f = 3/t_step_start
x = 1#np.sin(2*np.pi*f*t)
x = x*(step_amp + (1-step_amp)*(n>=t_step_start*Fs)*(n<=t_step_stop*Fs))

#=============================================================
#Basic limiter


#Define Limiter Threshold
thres = 0.5

#Estimate envelope
x_peak = gain_factor_smoothing(np.abs(x), 0.001, 0.04, Fs)

#gain function
g = np.minimum(np.ones(len(x_peak)), thres / x_peak)
x_lim = x * g

#=============================================================

fig, ax = plt.subplots(3, 1, sharex=True)

ax[0].plot(t, x, label=r'$x[n]$')
ax[0].plot(t, thres*np.ones_like(t), 'r--', label='Threshold')

ax[1].plot(t, np.abs(x), label=r'$|x[n]|$')
ax[1].plot(t, x_peak, label='Envelope')

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
    ax[i].set_xlim([0, t_max])

ax[2].legend(loc='lower left')

plt.show()
