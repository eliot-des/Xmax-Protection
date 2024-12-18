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
from filters import bilinear2ndOrder


#=============================================================
#Implement a limiter using a compensation filter based on physical approach of the speaker.

#Check limiter_approach_v1.png to see the signal's flowchart.
#=============================================================
#=============================================================

music = 'Thriller'
Fs, u = read(f'Audio/{music}.wav') 
u = u[:,0]          #select only one channel
u = normalize(u)    #normalize the signal to 1
G = 30              #gain of the amplifier -> Max tension in volts
u*=G                #tension in volts

t = np.arange(0, len(u)/Fs, 1/Fs)
tstart   = 0        #start time in seconds
duration = 4        #duration time in seconds

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

# bd_xu, ad_xu = sig.bilinear(b_xu, a_xu, fs=Fs)
bd_xu, ad_xu = bilinear2ndOrder(b_xu, a_xu, Fs)


#=============================================================

# Define Limiter Threshold
thres = Xmax_

# Envelope estimation parameters
attack_time   = 0.005
hold_time     = 0.005
release_time  = 0.010

release_coeff = 1 - np.exp(-2.2/(Fs*release_time))

# Convert times to samples
N_attack  = int(attack_time * Fs)
N_hold    = int(hold_time * Fs)


# Arrays to store results
u_hp    = np.zeros_like(u)      #tension supplying the speaker
x       = np.zeros_like(u)      #displacement
x_g     = np.ones_like(u)       #gain computer function
x_g_bis = np.ones_like(u)  
c       = np.ones_like(u)       #minimum filter and exp release for the computer gain
g       = np.ones_like(u)       #averaging
Cms_comp= np.ones_like(u)*Cms   #compensation compliance
gain_DC = np.zeros_like(u)      #gain @0Hz from compensation filter
gain_f  = np.zeros_like(u)      #broadband max gain of compensation filter



a_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])
b_comp = a_xu
#compute the digital coefficients
bd_comp, ad_comp = bilinear2ndOrder(b_comp, a_comp, Fs)


f = np.geomspace(10, 3000, 100)
Q0 = 1/np.sqrt(2)                        # neutral quality factor


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
    
    Cms_comp[n] = Cms*g[n]

    # Compare quality factor
    Q0_c = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms_comp[n])     # Quality factor for the compensation filter - c stands for compensation

    if Q0_c > Q0:
        Rms_comp = 1/Q0*np.sqrt(Mms/Cms_comp[n]) - (Bl**2/Rec)
    else:
        Rms_comp = Rms

    #we then compute the new compensation coeff filter (zeros are unchanged)...
    a_comp = np.array([Mms, Rms_comp+Bl**2/Rec, 1/Cms_comp[n]])

    #compute the digital coefficients
    bd_comp, ad_comp = bilinear2ndOrder(b_comp, a_comp, Fs)

    gain_DC[n] = 20*np.log10(np.sum(bd_comp)/np.sum(ad_comp))  #gives the gain in dB @0Hz

    # Search the max gain of the compensation filter
    _, H = sig.freqz(bd_comp, ad_comp, worN=f, fs=Fs)
    gain_f[n] = np.max(np.abs(H))


    #Apply the compensation filter to the delayed input tension signal
    u_hp[n] = bd_comp[0]*u[n-N_attack] + d_TDFII[0]

    d_TDFII[0] = bd_comp[1]*u[n-N_attack] - ad_comp[1]*u_hp[n] + d_TDFII[1]
    d_TDFII[1] = bd_comp[2]*u[n-N_attack] - ad_comp[2]*u_hp[n]



x_lim = sig.lfilter(bd_xu, ad_xu, u_hp)    #displacement in m

fig, ax = plt.subplots(3, 1, sharex=True)

#fig.suptitle(f'Sinusoidal signal of frequency {f:.2f} Hz with amplitude modulation.')

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, np.roll(u_hp, -N_attack), 'k--', label=r'$u_{hp}[n-N_{attack}]$')
ax[0].set(ylabel='Amplitude [V]')
ax[0].set_ylim(-G,G)

ax[1].plot(t, 20*np.log10(x_g), label='Gain computer output')
# ax[1].plot(t, x_g_bis, label='Gain computer bis output')
ax[1].plot(t, 20*np.log10(c), label=r'$c[n]$')
ax[1].plot(t, 20*np.log10(g), color='crimson', label=r'$g[n]$')
ax[1].plot(t, gain_DC ,'-.k', label=r'$gain_{DC}[n]$')
ax[1].plot(t, 20*np.log10(gain_f) ,'m', label=r'$max gain_{broadband}[n]$')
# ax[1].set_ylim([0, 1])


ax[2].plot(t, x, label=r'$x[n]$')
ax[2].plot(t, np.roll(x_lim, -N_attack),'k--', label=r'$x_{lim}[n-N_{attack}]$')
ax[2].plot(t, thres*np.ones_like(t), 'r--', label='+/-Threshold')
ax[2].plot(t, -thres*np.ones_like(t), 'r--')
ax[2].set(xlabel='Time [s]')

ax2twin = ax[2].twinx()
ax2twin.plot(t, g, 'g',alpha=0.2, label='Gain function')
ax2twin.set(ylabel='Gain')
ax2twin.legend(loc='upper right')

for i in range(3):
    ax[i].set(ylabel='Amplitude')
    ax[i].legend(loc='lower right')
    ax[i].grid()
    ax[i].set_xlim([t[0], t[-1]])
plt.show()

#=============================================================
# write tension signal u and u_hp to a wav file

# norm_factor = np.max(np.abs(u))

# u    = u   /norm_factor
# u_hp = u_hp/norm_factor

# u   *=32767
# u_hp*=32767

# write(f'Audio/Limiter/Approach_2/{music}_u.wav', Fs, u.astype(np.int16))
# write(f'Audio/Limiter/Approach_2/{music}_u_hp.wav', Fs, u_hp.astype(np.int16))