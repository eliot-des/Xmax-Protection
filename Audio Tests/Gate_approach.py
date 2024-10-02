import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize
import sys

Fs, v = read('Audio Tests/ShookOnesDisplacement.wav')
t = np.arange(0, len(v)/Fs, 1/Fs)
v = normalize(v)    #normalize the signal to 1

tstart = 1.9
tstop = 2.4

t = t[int(tstart*Fs):int(tstop*Fs)]
v = v[int(tstart*Fs):int(tstop*Fs)]**2

Xmax = 0.5
M = 0.2


# sys.exit()

Thigh = 0
T = 0
Tlow = 0
Thigh_return = 0
T_return = 0
Tlow_return = 0
Act = np.zeros_like(v)
Act_test = np.zeros_like(v)

for i in range(len(v)):
    if v[i] > Xmax-M:
        Tlow = 1
    if v[i] > Xmax and Tlow == 1:
        T = 1
    if v[i] > Xmax+M and T == 1 and Tlow == 1:
        Thigh = 1
        Tlow = 0
        T = 0
    if v[i] < Xmax+M and Thigh == 1:
        Thigh_return = 1
    if v[i] < Xmax and Thigh_return == 1:
        T_return = 1
        T = 0
    if v[i] < Xmax-M and T_return == 1 and T_return == 1:
        Tlow_return = 1
        Tlow = 0
        # Tlow = 0
        # Act_test[i] = 1
    if Thigh ==1 and Tlow == 0 and Tlow_return == 1:
        Thigh = 0
    Act[i]=Thigh


fig, ax = plt.subplots()
ax.plot(t,v, label='Squared signal')
ax.plot(t,Act, label='Activation func')
# ax.plot(t,Act_test, label='Activation func')
ax.axhline(y=Xmax+M, color='red', linestyle='--', label='Xmax treshold')
ax.axhline(y=Xmax, color='m', linestyle='-.', label='Xmax treshold -')
ax.axhline(y=Xmax-M, color='m', linestyle=':', label='Xmax treshold --')
ax.set_xlabel('Time [sec]')
# ax.set_ylim([0, 1+Xmax])
ax.legend()
ax.grid()
plt.show()