import numpy as np
import matplotlib.pyplot as plt

#=============================================================

#Implement a basic hard-clipping, for which the side chain is:

#-take the absolute value of the input signal
#-apply the gain computer function to the absolute value
#-apply a gain-smoothing function to the gain computer output

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

#=============================================================

# Define Limiter Threshold
thres = 0.5

# Envelope estimation parameters
attack_time   = 0.0001  # fast attack
release_time  = 0.04    # slower release
attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

# Convert times to samples
attack_samples = int(attack_time * Fs)
release_samples = int(release_time * Fs)

# Arrays to store results
x_lim = np.zeros_like(x)
x_g   = np.zeros_like(x)
g     = np.ones_like(x)


#=============================================================
# Implement the limiter in a for loop

for n in range(1, len(x)):

    # Calculate the absolute value of the input signal
    abs_x = np.abs(x[n])

    # Apply the gain computer function to the absolute value
    x_g[n] = np.minimum(1, thres/abs_x)

    # Apply the gain smoother to the gain computer output
    if x_g[n] < g[n-1]:
        k = attack_coeff
    else:
        k = release_coeff

    g[n] = (1 - k) * g[n-1] + k * x_g[n]

    # Apply the gain function
    x_lim[n] = x[n] * g[n]

#=============================================================

fig, ax = plt.subplots(3, 1, sharex=True)

ax[0].plot(t, x, label=r'$x[n]$')
ax[0].plot(t, thres*np.ones_like(t), 'r--', label='Threshold')


ax[1].plot(t, x_g, label='Gain computer output')
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
    ax[i].set_xlim([0, t_max])



plt.show()
