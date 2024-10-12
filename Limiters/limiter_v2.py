import numpy as np
import matplotlib.pyplot as plt

#=============================================================

#Implement a basic hard-clipping limiter, with a lookahead delay

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
x = np.sin(2*np.pi*f*t)
x = x*(step_amp + (1-step_amp)*(n>=t_step_start*Fs)*(n<=t_step_stop*Fs))

#=============================================================
#Basic limiter with lookahead delay


#Define Limiter Threshold
thres = 0.5

#Envelope estimation parameters
attack_time   = 0.001
release_time  = 0.04
attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

N = int((attack_time)*2.5* Fs) #delay in samples
print('delay in time:', N/Fs)

x_peak = np.zeros_like(x)
x_lim  = np.zeros_like(x)
g      = np.zeros_like(x)

for n in range(N, len(x)): #start from N to avoid problems with the buffer
    x_buff = x[n-N]

    #envelope estimation
    abs_x = np.abs(x[n])

    if abs_x > x_peak[n-1]:
        k = attack_coeff
    else:
        k = release_coeff

    x_peak[n] = (1 - k) * x_peak[n-1] + k * abs_x

    #gain function -> hard clipping
    g[n] = np.minimum(1, thres / x_peak[n])

    x_lim[n] = x_buff * g[n]


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



plt.show()
