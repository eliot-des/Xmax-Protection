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

#=============================================================
#Implement a basic hard-clipper.

#Check limiter_approach_v1.png to see the signal's flowchart.
#=============================================================
#=============================================================

music = 'Sacrifice2'
Fs, u = read(f'Audio/{music}.wav') 
u = u[:,0]          #select only one channel
u = normalize(u)    #normalize the signal to 1
G = 10              #gain of the amplifier -> Max tension in volts
u*=G             #tension in volts

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



# Low-frequency approximation of displacement
b_xu = np.array([0, 0, Bl/Rec])         
a_xu = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])

bdo_xu, ado_xu = sig.bilinear(b_xu, a_xu, fs=Fs)

# Just to save the filter coefficients before stabilization.
# o stands for original
bd_xu, ad_xu = bdo_xu, ado_xu 


#=============================================================
# Stabilization of the reciprocal filter, by putting the zeros
# of the X/U filter inside the unit circle

z, p, k = sig.tf2zpk(bd_xu, ad_xu)
#stabilization factor:
alpha = 0.9

for i in range(len(p)):
    if np.abs(z[i]) >= 1:
        z[i] = (1-alpha)/np.conjugate(z[i]) 

bd_xu, ad_xu = sig.zpk2tf(z, p, k)

f = np.geomspace(20, Fs/2, 1000)
_, H_analog = sig.freqs(b_xu,   a_xu, worN=f*2*np.pi)
_, H_digit  = sig.freqz(bd_xu, ad_xu, worN=f, fs=Fs)

# get module difference between analog and digital filter at the 
# lowest frequency, thencmultiply the gain of the digital filter 
# by delta_k to match the analog filter's modulus at low freqs.
delta_k = np.abs(H_analog[0])/np.abs(H_digit[0])
bd_xu, ad_xu = sig.zpk2tf(z, p, k*delta_k)

bd_ux = ad_xu
ad_ux = bd_xu

#normalize the reciprocal filter
ad0 = ad_ux[0]
ad_ux = ad_ux/ad0
bd_ux = bd_ux/ad0
#=============================================================

# Define Limiter Threshold
thres = Xmax_
knee  = 0.25

# Envelope estimation parameters
attack_time      = 0.004
hold_time        = 0.004
max_release_time = 1
average_time     = 0.2


average_coeff     = np.exp(-1/(Fs*average_time))
max_release_coeff = np.exp(-1/(Fs*max_release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
u_hp   = np.zeros_like(u)
x      = np.zeros_like(u)
x_lim  = np.zeros_like(x)
x_g    = np.ones_like(x)
c      = np.ones_like(x)
g      = np.ones_like(x)
u_peak = np.zeros_like(x)
u_rms  = np.zeros_like(x)
cf_u   = np.zeros_like(x)

d_xu = np.zeros(4) #memory for DF1 X/U filter 
d_ux = np.zeros(4) #memory for DF1 U/X filter
#=============================================================
# Implement the limiter in a for loop

for n in range(N_attack+N_hold, len(x)):

    #get crest factor of the tension signal
    
    abs_u = np.abs(u[n])
    u_peak[n] = np.maximum(abs_u, average_coeff*abs(u_peak[n-1]) + (1-average_coeff)*abs_u)
    u_rms[n]  = average_coeff*abs(u_rms[n-1]) + (1-average_coeff)*abs_u

    cf_u[n] = u_peak[n]/u_rms[n]
    release_time = 2*max_release_time/(cf_u[n]**2) - (2*attack_time/(cf_u[n]**2))
    release_coeff = np.exp(-1/(Fs*release_time))
    #release_coeff = max_release_coeff
    #filter the tension signal to get the displacement signal
    x[n] = bd_xu[0]*u[n-1] + bd_xu[1]*d_xu[0] + bd_xu[2]*d_xu[1] - ad_xu[1]*d_xu[2] - ad_xu[2]*d_xu[3]

    d_xu[1] = d_xu[0]
    d_xu[0] = u[n-1]
    d_xu[3] = d_xu[2]
    d_xu[2] = x[n]
    
    #x[n]*=1e3 #convert to mm
    # Calculate the absolute value of the input signal
    abs_x = np.abs(x[n])

    # Apply the gain computer function to the absolute value
    x_g[n] = np.minimum(1, thres/abs_x)
    #x_g[n] =  1/(1 + (abs_x/thres)**(1/knee))**knee
    
    # Apply the minimum filter to the gain computer output
    c[n] = np.min(x_g[n-(N_attack+N_hold):n+1])

    # Apply exponential release to the minimum filter output
    c[n] = np.minimum(c[n], release_coeff*c[n-1] +  (1-release_coeff)*c[n])

    # Apply the averaging filter to the minimum filter output
    g[n] = g[n-1] + 1/N_attack * (c[n] - c[n-N_attack])#recursive implementation

    # Apply the gain function to the delayed gain computer output
    x_lim[n] = x[n-N_attack] * g[n]
    
    #filter the limited displacement signal back to the tension signal
    u_hp[n]= bd_ux[0]*x_lim[n] + bd_ux[1]*d_ux[0] + bd_ux[2]*d_ux[1] - ad_ux[1]*d_ux[2] - ad_ux[2]*d_ux[3]

    d_ux[1] = d_ux[0]
    d_ux[0] = x_lim[n]
    d_ux[3] = d_ux[2]
    d_ux[2] = u_hp[n]


x_lim_true = sig.lfilter(bdo_xu, ado_xu, u_hp)

#=============================================================

fig, ax = plt.subplots(3, 1, sharex=True)

#fig.suptitle(f'Sinusoidal signal of frequency {f:.2f} Hz with amplitude modulation.')

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, np.roll(u_hp, -N_attack-1), 'k--', label=r'$u_{hp}[n]$')
ax[0].set(ylabel='Amplitude [V]')
ax[0].set_ylim(-G,G)

ax[1].plot(t, u_peak, label='Peak tension')
ax[1].plot(t, u_rms, label='RMS tension')
ax[1].plot(t, cf_u, 'r',label='Crest factor')
ax[1].set(ylim=[0,G], ylabel='Amplitude [V]')

ax[-1].plot(t, x, label=r'$x[n]$')
ax[-1].plot(t, np.roll(x_lim_true,-N_attack),'r', label=r'$x_{lim}^{(true)}[n+N_{attack}+1]$')
ax[-1].plot(t, np.roll(x_lim, -N_attack),'k--', label=r'$x_{lim}[n+N_{attack}]$')
ax[-1].plot(t, thres*np.ones_like(t), 'r--', label='Threshold')
ax[-1].plot(t,-thres*np.ones_like(t), 'r--')
ax[-1].set(xlabel='Time [s]')

ax2twin = ax[-1].twinx()
ax2twin.plot(t, g, 'g',alpha=0.2, label='Gain function')
ax2twin.set(ylabel='Gain')
ax2twin.legend(loc='upper right')

for i in range(3):
    ax[i].set(ylabel='Amplitude')
    ax[i].legend()
    ax[i].grid()
    ax[i].set_xlim([t[0], t[-1]])
plt.show()

#normalize u, u_hp by the normalization factor applied to u

'''
norm_factor = np.max(np.abs(u))

u=u/norm_factor
u_hp=u_hp/norm_factor

u*=32767
u_hp*=32767

#write(f'Audio/Limiter/Approach_1/{music}_u.wav', Fs, u.astype(np.int16))
write(f'Audio/Limiter/Approach_1/{music}_u_hp_AR.wav', Fs, u_hp.astype(np.int16))
'''