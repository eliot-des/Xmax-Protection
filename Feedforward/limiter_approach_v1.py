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

music = 'Thriller'
Fs, u = read(f'Audio/{music}.wav') 
u = u[:,0]          #select only one channel
u = normalize(u)    #normalize the signal to 1
G = 6.2               #gain of the amplifier -> Max tension in volts
u*=G             #tension in volts

t = np.arange(0, len(u)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 3            #duration time in seconds

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

bd_xu, ad_xu = sig.bilinear(b_xu, a_xu, fs=Fs)
bd_ux, ad_ux = sig.bilinear(a_xu, b_xu, fs=Fs)

'''
z, p, k = sig.tf2zpk(bd_ux, ad_ux)
# Add a small perturbation to the poles to ensure stability
p = p + 0.01 + 1j * 0.001
bd_ux, ad_ux = sig.zpk2tf(z, p, k)
'''
#=============================================================

# Define Limiter Threshold
thres = 1e-3

# Envelope estimation parameters
attack_time   = 0.002
hold_time     = 0.004
release_time  = 0.015

release_coeff = 1 - np.exp(-2.2/(Fs*release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
u_lim = np.zeros_like(u)
x     = np.zeros_like(u)
x_lim = np.zeros_like(x)
x_lim_= np.zeros_like(x)
x_g   = np.ones_like(x)
c     = np.ones_like(x)
g     = np.ones_like(x)

d_xu = np.zeros(4) #memory for DF1 X/U filter 
d_ux = np.zeros(4) #memory for DF1 U/X filter
#=============================================================
# Implement the limiter in a for loop
'''
x = sig.lfilter(bd_xu, ad_xu, u)    #displacement in m
x*= 1e3                             #displacement in mm
'''
for n in range(N_attack+N_hold, len(x)):
    
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
    #x_g[n] =  gain_computer_compressor(abs_x, thres, 1.5)

    # Apply the minimum filter to the gain computer output
    c[n] = np.min(x_g[n-(N_attack+N_hold):n+1])

    # Apply exponential release to the minimum filter output
    c[n] = np.minimum(c[n], (1-release_coeff)*c[n-1] + release_coeff*c[n]) #recursive implementation

    # Apply the averaging filter to the minimum filter output
    if n >= N_attack:
        g[n] = g[n-1] + 1/N_attack * (c[n] - c[n-N_attack]) #recursive implementation
    else:
        g[n] = c[n]  # No filtering for the initial samples

    # Apply the gain function to the delayed gain computer output
    x_lim[n] = x[n-N_attack] * g[n]
    
    #filter the limited displacement signal back to the tension signal
    u_lim[n]= bd_ux[0]*x_lim[n-1] + bd_ux[1]*d_ux[0] + bd_ux[2]*d_ux[1] - ad_ux[1]*d_ux[2] - ad_ux[2]*d_ux[3]

    d_ux[1] = d_ux[0]
    d_ux[0] = x_lim[n-1]
    d_ux[3] = d_ux[2]
    d_ux[2] = u_lim[n]

    #u_lim[n]/=1e3




u_lim2 = sig.lfilter(bd_ux, ad_ux, x_lim)


#=============================================================

fig, ax = plt.subplots(3, 1, sharex=True)

#fig.suptitle(f'Sinusoidal signal of frequency {f:.2f} Hz with amplitude modulation.')

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, np.roll(u_lim, -N_attack-2), 'k--', label=r'$u_{lim}[n]$')
ax[0].plot(t, np.roll(u_lim2, -N_attack-1), 'r--', label=r'$u_{lim2}[n]$')
ax[0].set(ylabel='Amplitude [V]')
ax[0].set_ylim(-G,G)
'''
ax[1].plot(t, x_g, label='Gain computer output')
ax[1].plot(t, c, label=r'$c[n]$')
ax[1].plot(t, g,color='crimson', label=r'$g[n]$')
'''

ax[2].plot(t, x, label=r'$x[n]$')
ax[2].plot(t, np.roll(x_lim, -N_attack),'k--', label=r'$x_{lim}[n]$')
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
    ax[i].set_xlim([t[0], t[-1]])
plt.show()


