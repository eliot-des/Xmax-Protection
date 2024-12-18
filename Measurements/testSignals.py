import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.io.wavfile import read, write

# Time
Fs = 48000
duration = 5
t = np.arange(0, duration, 1/Fs)

# Chirp
Fstart = 20
Fstop = 10000
chirp = sig.chirp(t, f0=Fstart, f1=Fstop, t1=duration, method='logarithmic', phi=90) # 
chirp = chirp[::-1]
write('Measurements/Chirp.wav', Fs, chirp.astype(np.float32))

# Sine
f0 = 50
Sine = np.sin(2*np.pi*f0*t)
# write('Measurements/Sine.wav', Fs, Sine.astype(np.float32))

# Plot
fig, ax = plt.subplots()
ax.plot(t, chirp)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Amp [a.u]')
plt.show()

# To mono a music
# signal_name = 'ShookOnes' 
# Fs, u = read('Audio/'+signal_name+'.wav')
# print("music dtype:", u.dtype)
# u = u[:,0]
# write('Measurements/'+signal_name+'_mono.wav', Fs, u.astype(np.int16))