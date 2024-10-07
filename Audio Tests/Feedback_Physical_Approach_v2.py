import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize, dynamic_peak_follower2, gain_factor_smoothing_sbs_bis, dynamic_peak_follower_sbs, gain_factor_smoothing_sbs
import sys

plt.rc('lines', linewidth=2)
plt.rc('font', size=14)
plt.rc('axes', linewidth=1.5, labelsize=14)
plt.rc('legend', fontsize=12)

#================================================================================
#============= Definition of a signal being the input voltage ===================
#================================================================================

# Fs, u = read('Audio Tests/Thriller.wav')
# t = np.arange(0, len(u[:,0])/Fs, 1/Fs)
# u = normalize(u[:,0])

#create a sweep/chirp signal
Fs = 48000
t = np.arange(0, 5, 1/Fs)
u = sig.chirp(t, f0=10000, f1=10, t1=5, method='logarithmic')

A = 10               #gain of the amplifier -> Max tension in volts
u = A*u             #tension in volts

tstart   = 0            #start time in seconds
duration = 5            #duration time in seconds

u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]

#================================================================================
#==== Importing speaker TS parameters and defining the displacement filter ======
#================================================================================

Xmax = 1        #in mm
Xmax_= Xmax/1e3 #in m
# Thiele-small parameters
loudspeakers = {'full-range1' :'Dayton_CE4895-8',
                'full-range2' :'Dayton_HARB252-8',
                'woofer1': 'Peerless_HDSP830860',
                'woofer2': 'Dayton_DCS165-4',
                'woofer3': 'Dayton_RS150-4',
                'subwoofer1': 'B&C_15FW76-4'}

loudspeaker = loudspeakers['full-range2']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)


#optimal minimal compliance to not 
R_v1 = (Xmax_*Rec/(A*Bl*Cms))

# Low-frequency approximation
B_LF = np.array([0, 0, Bl/Rec])         
A_LF = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])


# Analytical z-domain coeff from Low-frequency approximation
b[0] = B_LF[2]/Fs**2
b = np.zeros(3)                         
b[1] = 2*b[0]
b[2] = b[0]

a = np.zeros(3)
a[0] = A_LF[2]/Fs**2 + 2*A_LF[1]/Fs + 4*A_LF[0]
a[1] = 2*A_LF[2]/Fs**2 - 8*A_LF[0]
a[2] = A_LF[2]/Fs**2 + 4*A_LF[0] - 2*A_LF[1]/Fs

#normalize coeff
a0 = a[0]
a = a/a0
b = b/a0


#================================================================================
#========== Filter the input voltage to get an estimated displacement ===========
#================================================================================
x2 = np.zeros_like(u)        #displacement in m
x2 = sig.lfilter(b, a, u)    #displacement in m
x2 = x2*1000                #displacement in mm

#================================================================================
#====================== Compensation filter definition ==========================
#========================= Initializing parameters ==============================
#================================================================================

Cms_comp = np.zeros(len(u))    #compensation compliance to be adjusted sample by sample
Cms_max = Cms                  #initial compensation compliance equal (almost) to the real one of the speaker
Cms_min = 0.9*R*Cms              #minimum compensation compliance is 10% of the real Cms (CHOICE)
ratio = 0.6
Cms_comp[0] = Cms_max #ratio*Cms

B_comp = A_LF
A_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms_comp[0]])

b_comp = np.zeros(3)
a_comp = np.zeros(3)

b_comp[0] = A_LF[0]*4*Fs**2+A_LF[1]*2*Fs+A_LF[2]
b_comp[1] = -2*A_LF[0]*4*Fs**2+2*A_LF[2]
b_comp[2] = A_LF[0]*4*Fs**2-A_LF[1]*2*Fs+A_LF[2]

a_comp[0] = A_LF[0]*4*Fs**2+A_LF[1]*2*Fs+A_comp[2]
a_comp[1] = -2*A_LF[0]*4*Fs**2+2*A_comp[2]
a_comp[2] = A_LF[0]*4*Fs**2-A_LF[1]*2*Fs+A_comp[2]

#normalize coeff
a0_comp = a_comp[0]
a_comp = a_comp/a0_comp
b_comp = b_comp/a0_comp

f = np.geomspace(1, Fs/2, 1000)
w   = 2*np.pi*f

#================================================================================
#================== Main loop - feedback control algorithm ======================
#================================================================================

#Initializations

wr = 2*np.pi*fs
Q0 = 1/np.sqrt(2)                               # neutral quality factor
# Q0 = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms)     # original loudspeaker quality factor

x        = np.zeros(len(u))    #displacement in mm
x_lim    = np.zeros(len(u))    #limited displacement in mm
x_peak   = np.zeros(len(u))    #enveloppe estimator in mm

attack_peak    = 0.00005      # Attack time in seconds
release_peak   = 0.01         # Release time in seconds
attack_smooth  = 0.01         # Attack time for the gain smoothing function
release_smooth = 0.5         # Release time for the gain smoothing function


d_DFI = np.zeros(4)     # memory buffer for Direct Form I
d_TDFII = np.zeros(2)   # memory buffer for Transposed Direct Form I
a_comp = np.zeros(3)    # denumerator coefficient of the compensation filter in the z domain
u_hp = np.zeros_like(u)
track_Rms = np.zeros_like(u_hp)

for i in range(1, len(u_hp)):
    
    #we first estimate the displacement from the u_hp[i] sample
    x[i] = b[0]*u_hp[i-1] + b[1]*d_DFI[0] + b[2]*d_DFI[1] - a[1]*d_DFI[2] - a[2]*d_DFI[3]

    d_DFI[1] = d_DFI[0]
    d_DFI[0] = u_hp[i-1]
    d_DFI[3] = d_DFI[2]
    d_DFI[2] = x[i]

    x[i] *= 1000

    # Peak envelope estimation
    # x_peak[i] = dynamic_peak_follower_sbs(x[i], x_peak[i-1], attack_peak, release_peak, Fs) # x_peak in mm
    # Klippel envelope estimation
    # x_peak[i] = np.sqrt(x[i]**2 + (((x[i] - x[i-1]) * Fs)/wr)**2)

    #we then compare the peak of the displacement signal with Xmax
    # if np.abs(x[i]) > Xmax:
    if x_peak[i] > Xmax:
        Cms_target = Cms_min
    else:
        Cms_target = Cms_max

    # we then apply the gain smoothing function to the Cms_comp value of the compensation filter
    Cms_comp[i] = gain_factor_smoothing_sbs_bis(Cms_target, Cms_comp[i-1], attack_smooth, release_smooth, Fs)
    Rms_comp = 1/Q0*np.sqrt(Mms/Cms_comp[i]) - (Bl**2/Rec)
    track_Rms[i] = Rms_comp

    # #we then compute the new compensation coeff filter (zeros are unchanged)...
    A_comp = np.array([Mms, Rms_comp+Bl**2/Rec, 1/Cms_comp[i]])

    #convert to digital filter
    b_comp[0] = A_LF[0]*4*Fs**2+A_LF[1]*2*Fs+A_LF[2]
    b_comp[1] = -2*A_LF[0]*4*Fs**2+2*A_LF[2]
    b_comp[2] = A_LF[0]*4*Fs**2-A_LF[1]*2*Fs+A_LF[2]

    a_comp[0] = A_comp[0]*4*Fs**2+A_comp[1]*2*Fs+A_comp[2]
    a_comp[1] = -2*A_comp[0]*4*Fs**2+2*A_comp[2]
    a_comp[2] = A_comp[0]*4*Fs**2-A_comp[1]*2*Fs+A_comp[2]

    #normalize coeff
    a0_comp = a_comp[0]
    a_comp = a_comp/a0_comp
    b_comp = b_comp/a0_comp

    #... then apply it to the tension signal
    u_hp[i] = b_comp[0]*u[i] + d_TDFII[0]

    d_TDFII[0] = b_comp[1]*u[i] - a_comp[1]*u_hp[i] + d_TDFII[1]
    d_TDFII[1] = b_comp[2]*u[i] - a_comp[2]*u_hp[i]


#================================================================================
#================================== Main plot ===================================
#================================================================================

fig, ax = plt.subplots()

ax.plot(t, x, '-.b', linewidth=1, label=f'Estimated displacement (feedback)\n Attack {np.int16(attack_smooth*1000)} ms\n Release {np.int16(release_smooth*1000)} ms')
# ax.plot(t, x_peak, 'r-', linewidth=1, alpha=0.4, label='x_{peak}')
ax.plot(t, x2, linewidth=1, alpha=0.2, label='Displacement (no feedback)')
# ax.plot(t, u_hp, 'm', label='Limited displacement (feedback)')
ax.axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
ax.axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
ax.set_ylabel('Amplitude [mm]')
ax.set_xlabel('Time [sec]')
ax.legend(loc='lower left')
ax.grid()

ax1twin = ax.twinx()
# ax1twin.plot(t, Cms_comp*1000, 'k', linewidth=1, label='Cms compensation')
ax1twin.plot(t, track_Rms, 'k', linewidth=1, label='Rms compensation')
# ax1twin.set_ylim([Cms_min*990, Cms_max*1100])
ax1twin.legend(loc='upper left')
# ax1twin.set(ylabel='Compliance [mm/N]')
ax1twin.set(ylabel='Rms comp [kg/s]')

# plt.savefig(f'Figures/physical_approach/simulation_{loudspeaker}_settings3.pdf')
plt.show()