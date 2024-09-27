import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read


def dynamic_peak_follower1(x, attack_time, release_time, sample_rate):
    x_peak = 0.0
    x_peak_output = np.zeros_like(x)

    # Calculate attack and release coefficients based on provided times
    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * sample_rate))
    release_coeff = 1 - np.exp(-2.2 / (release_time * sample_rate))

    for n in range(len(x)):
        abs_x = np.abs(x[n])
        
        delta = abs_x - x_peak
        
        if delta > 0:
            x_peak += attack_coeff * delta
        else:
            x_peak += release_coeff * delta
        
        # Store the result in the output array
        x_peak_output[n] = x_peak

    return x_peak_output

def dynamic_peak_follower2(x, attack_time, release_time, Fs):
    
    x_peak = 0  
    x_peak_output = np.zeros_like(x)
    
    # Attack and release coefficients based on time constants
    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
    release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

    # Process each sample
    for n in range(len(x)):

        abs_x = np.abs(x[n])
        
        if abs_x > x_peak:
            x_peak = (1 - attack_coeff) * x_peak + attack_coeff * abs_x
        else:
            x_peak = (1 - release_coeff) * x_peak
        
                
        # Store the result
        x_peak_output[n] = x_peak
    
    return x_peak_output

#RMS Measurement
def dynamic_rms_follower(x, averaging_time, Fs):
    
    x2_rms = 0  
    x_rms_output = np.zeros_like(x)

    averaging_coeff = 1 - np.exp(-2.2 / (averaging_time * Fs))

    for n in range(len(x)):
        #squared RMS value
        x2_rms = (1 - averaging_coeff) * x2_rms + averaging_coeff * x[n]**2
        x_rms_output[n] = np.sqrt(x2_rms)
    
    return x_rms_output


def gain_factor_smoothing(x, attack_coeff, release_coeff):
    y = np.zeros_like(x)

    for n in range(1, len(x)):
        if x[n] > x[n-1]:
            k = attack_coeff
        else:
            k = release_coeff

        y[n] = (1 - k) * y[n-1] + k * x[n]

    return y


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
