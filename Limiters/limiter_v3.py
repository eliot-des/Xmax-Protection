import numpy as np
import matplotlib.pyplot as plt
#=============================================================

#Implement a basic hard-clipping limiter with:
#-a lookahead delay
#-a max filter to estimate the envelope

#Base on the paper:
# Smoothing of the Control Signal without Clipped Output in 
# Digital Peak Limiters
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

#Define Limiter Threshold
thres = 0.5


#Envelope estimation parameters
attack_time   = 0.0005
release_time  = 0.04

overshoot = 1.01
N = int((attack_time)*5* Fs) #delay in samples
beta = (1-overshoot)**(N+1)#1 - 1/overshoot
print(N)

a = 1 - 10**(np.log10((overshoot-1)/overshoot)/(N+1))


x_peak = np.zeros_like(x)
x_lim  = np.zeros_like(x)
e_max  = np.zeros_like(x)
x_max  = np.zeros_like(x)
g      = np.zeros_like(x)
c      = np.zeros_like(x)

for n in range(N, len(x)): #start from N to avoid problems with the buffer
    x_buff = x[n-N]

    #envelope estimation
    abs_x = np.abs(x[n])

    c[n] = np.maximum(abs_x, (abs_x - beta*e_max[n-1])/(1-beta))

    #apply statistic max filter on c from n-N to n
    x_max[n] = np.max(c[n-N:n+1])

    e_max[n] = (1 - a) * e_max[n-1] + a * x_max[n]

    #gain function -> hard clipping
    g[n] = np.minimum(1, thres / e_max[n])

    x_lim[n] = x_buff * g[n]


#=============================================================

fig, ax = plt.subplots(3, 1, sharex=True)

ax[0].plot(t, x, label=r'$x[n]$')
ax[0].plot(t, thres*np.ones_like(t), 'r--', label='Threshold')

ax[1].plot(t, np.abs(x), label=r'$|x[n]|$')
ax[1].plot(t, x_max, label='Envelope')

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
