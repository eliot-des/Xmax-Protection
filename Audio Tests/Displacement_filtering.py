import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.io.wavfile import read, write  

#import/read wavefile
Fs, data = read('Audio Tests/ShookOnes.wav')
data = data[:, 1]

# Thiele-small parameters
# Data from Peerless HDS-P830860 Datasheet

fs = 72
Rec = 6.4
Lec = 0.278e-3

Qms = 2.08
Qes = 0.725
Qts = 0.54

Mms = 8.88e-3
Cms = 560e-6
Rms = 1/(2*np.pi*fs*Cms*Qms)

Bl  = 5.74


#create filter
#analog coefficients
b = np.array([0, 0, 0, Bl])
a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])


#convert to digital filter
b, a = sig.bilinear(b, a, Fs)

#apply filter to signal

filtered_data = sig.lfilter(b, a, data)

fig, ax = plt.subplots(2, 1)

ax[0].plot(data)
ax[1].plot(filtered_data)

plt.show()


#normalize filtered data
filtered_data = filtered_data/np.max(np.abs(filtered_data))

filtered_data = filtered_data*2**15

write('Audio Tests/ShookOnesFiltered.wav', Fs, filtered_data.astype(np.int16))