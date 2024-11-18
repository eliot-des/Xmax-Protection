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
from filters import bilinear2ndOrder, EqualizerFilter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 20,
    "legend.fontsize": 18})

#=============================================================
#Implement a limiter using a LOW-SHELF filter.

#Check limiter_approach_v1.png to see the signal's flowchart.
#=============================================================
#=============================================================

music = 'Sacrifice1'
# music = 'ShookOnes_rearranged'
Fs, u = read(f'Audio/{music}.wav') 
u = u[:,0]          #select only one channel
u = normalize(u)    #normalize the signal to 1
G = 10               #gain of the amplifier -> Max tension in volts
u*=G             #tension in volts

t = np.arange(0, len(u)/Fs, 1/Fs)
tstart   = 0.5            #start time in seconds
duration = 2              #duration time in seconds

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

loudspeaker = loudspeakers['full-range1']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)



# Low-frequency approximation of displacement
b_xu = np.array([0, 0, Bl/Rec])         
a_xu = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])

bd_xu, ad_xu = sig.bilinear(b_xu, a_xu, fs=Fs)


#=============================================================

filter = EqualizerFilter("LS", fc=fs, Q=1/np.sqrt(2), dBgain=-6, fs=Fs)

# f = np.geomspace(1, Fs/2, 1000)
# w   = 2*np.pi*f

# _, H = sig.freqz(filter.b, filter.a, worN=f, fs=Fs)

# fig, ax = plt.subplots()
# ax.semilogx(f, 20*np.log10(np.abs(H)), 'b', label='Low-Shelf')
# ax.set_xlabel('Frequency  [Hz]')
# ax.set_ylabel('Gain [dB]')
# ax.grid(which='both')
# ax.legend(loc='lower right')
# plt.show()
# sys.exit()

#=============================================================

# Define Limiter Threshold
thres = Xmax_

# Envelope estimation parameters
attack_time   = 0.003
hold_time     = 0.010
release_time  = 0.085

release_coeff = 1 - np.exp(-2.2/(Fs*release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
u_hp    = np.zeros_like(u)      #tension supplying the speaker
x       = np.zeros_like(u)      #displacement
x_g     = np.ones_like(u)       #gain computer function
x_g_bis = np.ones_like(u)  
c       = np.ones_like(u)       #gain computer after the minimum filter output
g       = np.ones_like(u)       #averaging
dBgain  = np.ones_like(u)       #LS filter gain



d_xu = np.zeros(4)   #memory for DF1 X/U filter 
d_TDFII = np.zeros(2) #memory for TDF2 U/X filter
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
    abs_u = np.abs(u[n])

    # Apply the gain computer function to get the CMS ratio allowing to not exceed Xmax
    # x_g[n] = np.minimum(1, thres*Rec/(abs_u*Bl*Cms))
    x_g[n] = np.minimum(1, thres/(np.abs(x[n])))
    
    # Apply the minimum filter to the gain computer output
    c[n] = np.min(x_g[n-(N_attack+N_hold):n+1])

    # Apply exponential release to the minimum filter output
    c[n] = np.minimum(c[n], (1-release_coeff)*c[n-1] + release_coeff*c[n])

    # Apply the averaging filter to the minimum filter output
    if n >= N_attack:
        g[n] = g[n-1] + 1/N_attack * (c[n] - c[n-N_attack]) #recursive implementation
    else:
        g[n] = c[n]  # No filtering for the initial samples
    
    dBgain[n] = 20*np.log10(g[n])

    #compute the digital coefficients of LS filter
    filter = EqualizerFilter("LS", fc=200, Q=1/np.sqrt(2), dBgain=dBgain[n], fs=Fs)

    #Apply the compensation filter to the delayed input tension signal
    u_hp[n] = filter.b[0]*u[n-N_attack] + d_TDFII[0]

    d_TDFII[0] = filter.b[1]*u[n-N_attack] - filter.a[1]*u_hp[n] + d_TDFII[1]
    d_TDFII[1] = filter.b[2]*u[n-N_attack] - filter.a[2]*u_hp[n]



x_lim = sig.lfilter(bd_xu, ad_xu, u_hp)    #displacement in m

fig, ax = plt.subplots(2, 1, figsize=(12,7), sharex=True, layout='constrained')

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, np.roll(u_hp, -N_attack-1), 'k--', label=r'$u_{lim}[n+N_{attack}]$')
ax[0].set(ylabel='Amplitude [V]')
ax[0].set_ylim(-G,G)
ax[0].legend(loc='lower left')

ax2twin = ax[-1].twinx()
ax2twin.plot(t, g, color='crimson',alpha=0.2, label='Gain function')
ax2twin.set(ylabel='Gain')
ax2twin.legend(loc='upper right')

ax[-1].plot(t, x*1e3, label=r'$x[n]$')
#ax[-1].plot(t, np.roll(x_lim, -N_attack)*1e3,'k', label=r'$x_{lim}[n+N_{attack}]$')
ax[-1].plot(t, np.roll(x_lim,-N_attack)*1e3,'k', label=r'$x_{lim}[n+N_{attack}]$') #label=r'$x_{lim}^{(true)}[n+N_{attack}]$'
ax[-1].plot(t, thres*np.ones_like(t)*1e3, 'r--', label='Threshold')
ax[-1].plot(t, -thres*np.ones_like(t)*1e3, 'r--')
ax[-1].set(xlabel='Time [s]',ylabel='Amplitude [mm]')
ax[-1].legend()


for i in range(len(ax)):
    ax[i].grid()
    ax[i].set_xlim([t[0], t[-1]])
plt.show()

fig.savefig("LowShelf_Sacrice.pdf", format="pdf", transparent=True)
