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
#=============================================================

music = 'Sacrifice2'
Fs, u = read(f'Audio/{music}.wav') 
u = u[:,0]          #select only one channel
u = normalize(u)    #normalize the signal to 1
G = 10               #gain of the amplifier -> Max tension in volts
u*=G             #tension in volts

t = np.arange(0, len(u)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 4            #duration time in seconds

u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]

'''
t_max = 0.1
t_step_start = t_max/3
t_step_stop  = 2*t_max/3
step_amp = 0.1

Fs = 44100
t = np.arange(0, t_max, 1/Fs)
n = np.arange(0, len(t))

f = 2/t_step_start
u = G*(step_amp + (1-step_amp)*(n>=t_step_start*Fs)*(n<=t_step_stop*Fs))
'''
#=============================================================

Xmax = 1.5          #in mm
Xmax_= Xmax*1e-3    #in m

# Thiele-small parameters
loudspeakers = {'full-range1' :'Dayton_CE4895-8',
                'full-range2' :'Dayton_HARB252-8',
                'woofer1': 'Peerless_HDSP830860',
                'woofer2': 'Dayton_DCS165-4',
                'woofer3': 'Dayton_RS150-4',
                'subwoofer1': 'B&C_15FW76-4'}

loudspeaker = loudspeakers['full-range1']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)



# Low-frequency approximation of X/U(s)
b_xu = np.array([0, 0, Bl/Rec])         
a_xu = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])

bd_xu, ad_xu = sig.bilinear(b_xu, a_xu, fs=Fs)

#=============================================================

# Define Limiter Threshold
thresh = Xmax_*0.8
knee = 0.1

# Envelope estimation parameters
attack_time   = 0.004
hold_time     = 0.008
release_time  = 0.015

release_coeff = 1 - np.exp(-2.2/(Fs*release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
u_hp    = np.zeros_like(u)
x       = np.zeros_like(u)
x_lim   = np.zeros_like(x)
x_g     = np.ones_like(u)
c       = np.ones_like(u)
g       = np.ones_like(u)


d_xu = np.zeros(4)   #memory for DF1 X/U filter
#=============================================================
# Implementation of the algorithm

for n in range(N_attack+N_hold, len(x)):
    #filter the tension signal to get the displacement signal
    x[n] = bd_xu[0]*u[n-1] + bd_xu[1]*d_xu[0] + bd_xu[2]*d_xu[1] - ad_xu[1]*d_xu[2] - ad_xu[2]*d_xu[3]

    d_xu[1] = d_xu[0]
    d_xu[0] = u[n-1]
    d_xu[3] = d_xu[2]
    d_xu[2] = x[n]
    
    # Calculate the absolute value of the input signal
    abs_x = np.abs(x[n])

    # Apply the gain computer function to the absolute value of the displacement
    #x_g[n] = np.minimum(1, thresh/abs_x)
    x_g[n] = 1/(1 + (abs_x/thresh)**(1/knee))**knee

    # Apply the minimum filter to the gain computer output
    c[n] = np.min(x_g[n-(N_attack+N_hold):n+1])

    # Apply exponential release to the minimum filter output
    c[n] = np.minimum(c[n], (1-release_coeff)*c[n-1] + release_coeff*c[n])

    # Apply the averaging filter to the minimum filter output

    #g[n] = np.dot(fir_coeffs, c[n-N_attack+1:n+1])
    g[n] = g[n-1] + 1/N_attack * (c[n] - c[n-N_attack]) #recursive implementation


    # Apply the gain function to the delayed tension and displacement signals to see what happends...
    x_lim[n] = x[n-N_attack] * g[n]
    u_hp[n]  = u[n-N_attack] * g[n]



#displacement signal estimated from the tension signal going to the loudspeaker
x_lim2 = sig.lfilter(bd_xu, ad_xu, u_hp)    


fig, ax = plt.subplots(3, 1, sharex=True)

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, np.roll(u_hp, -N_attack), 'k--', label=r'$u_{hp}[n+N_{attack}]$')
ax[0].set(ylabel='Amplitude [V]')
ax[0].set_ylim(-G,G)

ax[1].plot(t, x_g, label='Gain computer output')
ax[1].plot(t, c, label=r'$c[n]$')
ax[1].plot(t, np.roll(g, -N_attack),color='crimson', label=r'$g[n+N_{attack}]$')

ax[2].plot(t, x, label=r'$x[n]$')
ax[2].plot(t, np.roll(x_lim2,-N_attack),'r', label=r'$x_{lim2}[n+N_{attack}]$')
ax[2].plot(t, np.roll(x_lim,-N_attack),'k--', label=r'$x_{lim}[n+N_{attack}]$')

ax[2].plot(t, thresh*np.ones_like(t), 'k--',alpha=0.5, label='Threshold')
ax[2].plot(t, -thresh*np.ones_like(t), 'k--',alpha=0.5)
ax[2].set(xlabel='Time [s]')

ax2twin = ax[2].twinx()
ax2twin.plot(t, np.roll(g, -N_attack), 'g',alpha=0.2, label='g[n+N_{attack}]')
ax2twin.set(ylabel='Gain')
ax2twin.legend(loc='upper right')

for i in range(3):
    ax[i].set(ylabel='Amplitude')
    ax[i].legend(loc='lower left')
    ax[i].grid()
    ax[i].set_xlim([t[0], t[-1]])
plt.show()

#=============================================================
# write tension signal u and u_hp to a wav file

'''
norm_factor = np.max(np.abs(u))

u    = u   /norm_factor
u_hp = u_hp/norm_factor

u   *=32767
u_hp*=32767

write(f'Audio/Limiter/Approach_3/{music}_u.wav', Fs, u.astype(np.int16))
write(f'Audio/Limiter/Approach_3/{music}_u_hp.wav', Fs, u_hp.astype(np.int16))
'''


fig, ax = plt.subplots(1, 2, sharey=True, layout='constrained')
ax[0].plot(x, u, 'o', alpha=0.1)
ax[0].set(xlabel='Displacement [m]', ylabel='Voltage [V]')
ax[1].plot(x_lim2, u_hp, 'o', alpha=0.1)

ax[1].set(xlabel='Displacement [m]')
#get ylimits of axe ax[1]:
ymin, ymax = ax[1].get_ylim()
ax[1].plot([-Xmax_, -Xmax_], [ymin, ymax], 'r')
ax[1].plot([Xmax_, Xmax_], [ymin, ymax], 'r')
ax[1].plot([-thresh, -thresh], [ymin, ymax], 'k--')
ax[1].plot([thresh, thresh], [ymin, ymax], 'k--')
ax[0].grid()
ax[1].grid()
plt.show()