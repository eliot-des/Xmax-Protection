import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.io.wavfile import read, write
import lib_impulse as lib

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
Fs, u = read(f'Audio/Thriller.wav') 
u = normalize(u)    #normalize the signal to 1
u = u[:,0]          #select only the left channel
G = 1               #gain of the amplifier -> Max tension in volts
u*= G               #tension in volts

t = np.arange(0, len(u)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 4            #duration time in seconds

u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
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

# Low-frequency approximation of displacement from tension
b_xu = np.array([0, 0, Bl/Rec])  
a_xu = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])

bd_xu, ad_xu = sig.bilinear(b_xu, a_xu, fs=Fs)

#=============================================================
#test to stabilize the future inverse filter of the displacement/tension filter
z, p, k = sig.tf2zpk(bd_xu, ad_xu)

#stabilization factor:
alpha = 0.9

for i in range(len(p)):
    if np.abs(z[i]) >= 1:
        # homemade formula, don't know is there exists something like this in the litterature
        z[i] = (1-alpha)/np.conjugate(z[i]) 
        print(f'z[{i}] = {z[i]}')

bd_xu, ad_xu = sig.zpk2tf(z, p, k)

_, H_analog = sig.freqs(b_xu,   a_xu, worN=f*2*np.pi)
_, H_digit  = sig.freqz(bd_xu, ad_xu, worN=f, fs=Fs)

# get module difference between analog and digital filter at the lowest frequency, then
# multiply the gain of the digital filter by delta_k to match the analog filter's modulus at low freqs.

delta_k = np.abs(H_analog[0])/np.abs(H_digit[0])
bd_xu, ad_xu = sig.zpk2tf(z, p, k*delta_k)
#=============================================================


_, H_digit  = sig.freqz(bd_xu, ad_xu, worN=f, fs=Fs)

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
# inverse filter to get tension from displacement, with the stabilization factor

# make the reciprocal filter to get tension from displacement
bd_ux = ad_xu
ad_ux = bd_xu

#normalize the reciprocal filter
ad0 = ad_ux[0]
ad_ux = ad_ux/ad0
bd_ux = bd_ux/ad0

# Process an audio signal, as being the tension signal, 
# then deduce the displacement signal, and then the tension signal 
# again to see the effect of the stabilization factor

x = sig.lfilter(bd_xu, ad_xu, u)
ub = sig.lfilter(bd_ux, ad_ux, x)



fig, ax = plt.subplots(3, 1, sharex=True)

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, ub,'--', label=r'$u_{b}[n]$')
ax[1].plot(t, x, label=r'$x_[n]$')
ax[2].plot(t, 20*np.log10(np.abs(u - ub)), label=r'$|u[n] - u_{b}[n]|$')


ax[0].set(ylabel='Amplitude [V]')
ax[1].set(ylabel='Amplitude [m]')
ax[2].set(xlabel='Time [s]', ylabel='Error [dB ref 1]')

for i in range(3):
    ax[i].legend()
    ax[i].grid()

plt.show()





