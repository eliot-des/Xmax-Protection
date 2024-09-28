import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read
from lib import dynamic_peak_follower1, dynamic_peak_follower2, dynamic_rms_follower


Fs = 48000
t  = np.arange(0, 0.19, 1/Fs)

x = np.zeros(len(t))
t_start = 0.02
t_stop  = 0.04

x[int(t_start*Fs): int(t_stop*Fs)] = 1

# Set parameters for the dynamic follower
attack_time    = 0.01    # Attack time in seconds
release_time   = 0.08    # Release time in seconds

x_peak1 = dynamic_peak_follower1(x, attack_time, release_time, Fs)
x_peak2 = dynamic_peak_follower2(x, attack_time, release_time, Fs)

# Plot the original and the peak-controlled signal
fig, ax = plt.subplots()

ax.plot(t, x, label=r'$x$')
ax.plot(t, x_peak1, label=r'$x_{peak1}$')
ax.plot(t, x_peak2, label=r'$x_{peak2}$')
ax.set(xlabel='Time [s]', ylabel='Amplitude', title='Peak-controlled signal')
ax.legend()
ax.grid()



Fs, x = read('Audio Tests/ShookOnes.wav')
x = x[:,0]
x = x/np.max(np.abs(x))
t = np.arange(len(x))/Fs

# Set parameters for the dynamic follower
attack_time    = 0.00005    # Attack time in seconds
release_time   = 0.07    # Release time in seconds
averaging_time = 0.02    # Averaging time in seconds

x_peak = dynamic_peak_follower1(x, attack_time, release_time, Fs)
x_rms  = dynamic_rms_follower(x, averaging_time, Fs)

# Plot the original and the peak-controlled signal
fig, ax = plt.subplots()
ax.plot(t, x, label=r'$x$')
ax.plot(t, x_peak,'k-', label=r'$x_{peak}$')
ax.plot(t, -x_peak,'k-')
ax.plot(t, x_rms, label=r'$x_{rms}$')
ax.set(xlabel='Time [s]', ylabel='Amplitude', title='Peak and rms follower signal')
ax.legend()
ax.grid()
plt.show()
