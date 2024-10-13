import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.io.wavfile import read, write

import os
import sys
# insert root directory into python module search path
fpath = os.path.join(os.getcwd(), 'Modules')
sys.path.append(fpath)

from audio import normalize

def bilinear2ndOrder(b, a, Fs):
    bd = np.zeros(3) 
    bd[0] = b[0]*4*Fs**2 + b[1]*2*Fs + b[2]
    bd[1] = -2*b[0]*4*Fs**2 + 2*b[2]
    bd[2] = b[0]*4*Fs**2 - b[1]*2*Fs + b[2]

    ad = np.zeros(3)
    ad[0] = a[0]*4*Fs**2 + a[1]*2*Fs + a[2]
    ad[1] = -2*a[0]*4*Fs**2 + 2*a[2]
    ad[2] = a[0]*4*Fs**2 - a[1]*2*Fs + a[2]

    #normalize coeff
    ad0 = ad[0]
    ad = ad/ad0
    bd = bd/ad0

    return bd, ad


#=============================================================
#Implement a basic hard-clipper.

#Check limiter_approach_v1.png to see the signal's flowchart.
#=============================================================
#=============================================================

music = 'Thriller'
Fs, x = read(f'Audio/Limiter/{music}.wav') 
x = normalize(x)    #normalize the signal to 1
G = 10               #gain of the amplifier -> Max tension in volts
x*=G             #tension in volts

t = np.arange(0, len(x)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 3            #duration time in seconds

x = x[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]

#=============================================================

Xmax = 1            #in mm
Xmax_= Xmax*1e-3    #in m

# Thiele-small parameters
loudspeakers = {'full-range1' :'Dayton_CE4895-8',
                'full-range2' :'Dayton_HARB252-8',
                'woofer1': 'Peerless_HDSP830860',
                'woofer2': 'Dayton_DCS165-4',
                'woofer3': 'Dayton_RS150-4',
                'subwoofer1': 'B&C_15FW76-4'}

loudspeaker = loudspeakers['woofer1']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)

f = np.geomspace(20, Fs/2, 1000)

# Low-frequency approximation of tension
b_ux = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])
a_ux = np.array([0, 0, Bl/Rec])  

b_ux = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])
a_ux = np.array([0, 0, 0, Bl])


bd_ux, ad_ux = sig.bilinear(b_ux, a_ux, fs=Fs)

z, p, k = sig.tf2zpk(bd_ux, ad_ux)
print(p)
p = p + 0.02
z = np.append(z, -1)
print(p)
bd_ux, ad_ux = sig.zpk2tf(z, p, k)



_, H_analog = sig.freqs(b_ux, a_ux, worN=f*2*np.pi)
_, H_digit = sig.freqz(bd_ux, ad_ux, worN=f, fs=Fs)





fig, ax = plt.subplots(2, 1,sharex=True)
ax[0].semilogx(f, 20*np.log10(np.abs(H_digit)), label='Digital')
ax[0].semilogx(f, 20*np.log10(np.abs(H_analog)), label='Analog')
ax[1].semilogx(f, np.angle(H_digit), label='Digital')
ax[1].semilogx(f, np.angle(H_analog), label='Analog')

ax[0].set(xlim=[20, 20000])
for i in range(2):
    ax[i].grid(which='both', axis='both')
    ax[i].legend()
plt.show()  


#=============================================================
#print tension signal

u = sig.lfilter(bd_ux, ad_ux, x)

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(t, x)
ax[1].plot(t, u)

plt.show()


