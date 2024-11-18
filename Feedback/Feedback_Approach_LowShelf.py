import os
import sys
# insert root directory into python module search path
fpath = os.path.join(os.getcwd(), 'Modules')
sys.path.append(fpath)

from audio import normalize
from filters import bilinear2ndOrder, EqualizerFilter


#=============================================================
#Implementation of an hybrid feedback/feedforward lowshelf filter
# to limit displacement of a loudspeaker.
#=============================================================
#=============================================================

music = 'ShookOnes'
# music = 'ShookOnes_rearranged'
Fs, u = read(f'Audio/{music}.wav') 
u = u[:,0]          #select only one channel
u = normalize(u)    #normalize the signal to 1
G = 1000               #gain of the amplifier -> Max tension in volts
u*=G             #tension in volts

t = np.arange(0, len(u)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 4              #duration time in seconds

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
u_hp_forward = np.zeros_like(u)      #tension supplying to the speaker throught the feedforward path
u_hp_feedback= np.zeros_like(u)      #tension supplying to the speaker obtained from the feedback path
x       = np.zeros_like(u)      #displacement
x_g     = np.ones_like(u)       #gain computer function
x_g_bis = np.ones_like(u)  
c       = np.ones_like(u)       #gain computer after the minimum filter output
g       = np.ones_like(u)       #averaging
dBgain  = np.ones_like(u)       #LS filter gain



d_xu = np.zeros(4)   #memory for DF1 X/U filter 
d_TDFII = np.zeros(2) #memory for TDF2 U/X filter
#=============================================================
for n in range(N_attack+N_hold, len(x)):

    #compute the digital coefficients of LS filter
    filter = EqualizerFilter("LS", fc=200, Q=1/np.sqrt(2), dBgain=dBgain[n], fs=Fs)

    #Apply the compensation filter to the delayed input tension signal
    u_hp_feedback[n] = filter.b[0]*u[n-N_attack] + d_TDFII[0]

    d_TDFII[0] = filter.b[1]*u[n-N_attack] - filter.a[1]*u_hp_feedback[n] + d_TDFII[1]
    d_TDFII[1] = filter.b[2]*u[n-N_attack] - filter.a[2]*u_hp_feedback[n]


    #filter the tension signal to get the displacement signal
    x[n] = bd_xu[0]*u_hp_feedback[n-1] + bd_xu[1]*d_xu[0] + bd_xu[2]*d_xu[1] - ad_xu[1]*d_xu[2] - ad_xu[2]*d_xu[3]

    d_xu[1] = d_xu[0]
    d_xu[0] = u_hp_feedback[n-1]
    d_xu[3] = d_xu[2]
    d_xu[2] = x[n]
    
    #x[n]*=1e3 #convert to mm
    # Calculate the absolute value of the input signal

    # Apply the gain computer function to get the CMS ratio allowing to not exceed Xmax
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

    



x_lim = sig.lfilter(bd_xu, ad_xu, u_hp)    #displacement in m

fig, ax = plt.subplots(3, 1, sharex=True)

ax[0].plot(t, u, label=r'$u[n]$')
ax[0].plot(t, np.roll(u_hp, -N_attack), 'k--', label=r'$u_{hp}[n-N_{attack}]$')
ax[0].set(ylabel='Amplitude [V]')
ax[0].set_ylim(-G,G)

ax[1].plot(t, x_g, label='Gain computer output')
ax[1].plot(t, c, label=r'$c[n]$')
ax[1].plot(t, g, color='crimson', label=r'$g[n]$')
ax[1].set_ylim([0, 1])


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
    ax[i].legend()
    ax[i].grid()
    ax[i].set_xlim([t[0], t[-1]])
plt.show()

'''
norm_factor = np.max(np.abs(u))

u=u/norm_factor
u_hp=u_hp/norm_factor

u*=32767
u_hp*=32767

write(f'Audio/FeedBack/Approach_LowShelf/{music}_u.wav', Fs, u.astype(np.int16))
write(f'Audio/FeedBack/Approach_LowShelf/{music}_u_hp.wav', Fs, u_hp.astype(np.int16))
'''