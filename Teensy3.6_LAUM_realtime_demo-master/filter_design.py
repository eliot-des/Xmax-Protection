import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

fs = 48000

fc = 200  # cut off frequency
b, a = signal.butter(2, fc/(fs/2), 'high')

f_axis = np.logspace(np.log10(20), np.log10(20e3), 1000)
_, H = signal.freqz(b, a, worN=f_axis, fs=fs)

fig, ax = plt.subplots()
ax.semilogx(f_axis, 20*np.log10(np.abs(H)))
ax.set(xlim=(20, 20e3))

print(f"volatile float b0 = {b[0]}, b1 = {b[1]}, b2 = {b[2]};")
print(f"volatile float a1 = {a[1]}, a2 = {a[2]};")

plt.show()
