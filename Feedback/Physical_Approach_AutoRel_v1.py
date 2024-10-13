import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
import os
import sys
# insert root directory into python module search path
fpath = os.path.join(os.getcwd(), 'modules')
sys.path.append(fpath)
from audio import normalize
from dynamic_trackers import dynamic_peak_follower2, gain_factor_smoothing_sbs_bis, dynamic_peak_follower_sbs, gain_factor_smoothing_sbs, rms_follower_sbs

plt.rc('lines', linewidth=2)
plt.rc('font', size=14)
plt.rc('axes', linewidth=1.5, labelsize=14)
plt.rc('legend', fontsize=12)


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

#================================================================================
#============= Definition of a signal being the input voltage ===================
#================================================================================

Fs, u = read('Audio/ShookOnes_rearrangedBis_Displacement.wav')
t = np.arange(0, len(u[:,0])/Fs, 1/Fs)
u = normalize(u[:,0])

# #create a sweep/chirp signal
# Fs = 48000
# t = np.arange(0, 5, 1/Fs)
# u = sig.chirp(t, f0=10000, f1=10, t1=5, method='logarithmic')

A = 4               #gain of the amplifier -> Max tension in volts
u = A*u             #tension in volts

# tstart   = 0            #start time in seconds
# duration = 5            #duration time in seconds
# tstart   = 12.14            #start time in seconds
# duration = 2.5            #duration time in seconds

# u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
# t = t[int(tstart*Fs):int((tstart+duration)*Fs)]

# fig, ax = plt.subplots()
# ax.plot(t, u, label='Signal')
# ax.set_xlabel('Time [sec]')
# ax.set_ylabel('Amplitude [u.a]')
# ax.grid()
# plt.show()

# sys.exit()

#================================================================================
#==== Importing speaker TS parameters and defining the displacement filter ======
#================================================================================

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
B_LF = np.array([0, 0, Bl/Rec])         
A_LF = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])


# Analytical z-domain coeff from Low-frequency approximation

b, a = bilinear2ndOrder(B_LF, A_LF, Fs)

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

R = Xmax_*Rec/(A*Bl*Cms)
print('Chosen R = ', np.round(R, 3))

Cms_comp = np.zeros(len(u))         #compensation compliance to be adjusted sample by sample
Cms_comp[0] = Cms                   #initial compliance is the one of the speaker
Cms_min = 0.9*R*Cms                 #minimum compensation with a real ratio of 110% of R (CHOICE)


B_comp = A_LF
A_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms_comp[0]])

b_comp, a_comp = bilinear2ndOrder(B_comp, A_comp, Fs)

f = np.geomspace(1, Fs/2, 1000)
w   = 2*np.pi*f

#================================================================================
#================== Main loop - feedback control algorithm ======================
#================================================================================

#Initializations

wr = 2*np.pi*fs                                 # resonance frequency of the speaker
Q0 = 1/np.sqrt(2)                               # neutral quality factor
# Q0 = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms)     # original loudspeaker quality factor
T_avg = 0.5                                     # average time for RMS estimator
Threshold = 0.07

x        = np.zeros(len(u))    #displacement in mm
x_peak   = np.zeros(len(u))    #enveloppe estimator in mm
rms      = np.zeros(len(u))    #RMS follower

attack_peak    = 0.00005      # Attack time in seconds
release_peak   = 0.01         # Release time in seconds
attack_smooth  = 0.01         # Attack time for the gain smoothing function
release_smooth = 0.5          # Release time for the gain smoothing function


d_DFI   = np.zeros(4)     # memory buffer for Direct Form I
d_TDFII = np.zeros(2)   # memory buffer for Transposed Direct Form I
u_hp    = np.zeros_like(u)
track_Rms = np.zeros_like(u_hp)
track_AutoR = np.zeros_like(u_hp)

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

    # RMS estimator
    rms[i] = rms_follower_sbs(x[i], rms[i-1], T_avg, Fs)

    if rms[i] > Threshold:
        AutoR = 50
    else:
        AutoR = 1
    
    track_AutoR[i] = AutoR

    # Klippel envelope estimation
    dx = (x[i] - x[i-1]) * Fs
    x_peak[i] = np.sqrt(x[i]**2 + (dx/wr)**2)

    #we then compare the peak of the displacement signal with Xmax
    if np.abs(x[i]) > Xmax:
    # if x_peak[i] > Xmax:
        Cms_target = Cms_min
    else:
        Cms_target = Cms

    # we then apply the gain smoothing function to the Cms_comp value of the compensation filter
    Cms_comp[i] = gain_factor_smoothing_sbs_bis(Cms_target, Cms_comp[i-1], attack_smooth, AutoR*release_smooth, Fs)
    Rms_comp = 1/Q0*np.sqrt(Mms/Cms_comp[i]) - (Bl**2/Rec)      # to impose the Q0 defined before
    Rms_comp = Rms                                              # Not the adjust Q0 (test)
    track_Rms[i] = Rms_comp

    #we then compute the new compensation coeff filter (zeros are unchanged)
    A_comp = np.array([Mms, Rms_comp+Bl**2/Rec, 1/Cms_comp[i]])

    #compute the digital coefficients
    b_comp, a_comp = bilinear2ndOrder(B_comp, A_comp, Fs)
    

    #Apply the compensation filter to the input tension signal
    u_hp[i] = b_comp[0]*u[i] + d_TDFII[0]

    d_TDFII[0] = b_comp[1]*u[i] - a_comp[1]*u_hp[i] + d_TDFII[1]
    d_TDFII[1] = b_comp[2]*u[i] - a_comp[2]*u_hp[i]


#================================================================================
#================================== Main plot ===================================
#================================================================================

fig, ax = plt.subplots(2,1, sharex=True)

ax[0].plot(t, x, 'b', linewidth=1, label=f'Estimated displacement (feedback)\n Attack {attack_smooth*1e3:.0f} ms\n Release {release_smooth*1e3:.0f} ms')
# ax[0].plot(t, x_peak, 'r-', linewidth=1, alpha=0.4, label='x_{peak}')
ax[0].plot(t, x2, linewidth=1, alpha=0.2, label='Displacement (no feedback)')
# ax[0].plot(t, u_hp, 'm', label='Limited displacement (feedback)')
ax[0].plot(t, rms, 'm', linewidth=2, alpha=0.6, label=f'RMS estimator (no feedback)\n Tavg {np.int16(T_avg*1000)} ms')
ax[0].axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
ax[0].axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
ax[0].set_ylabel('Amplitude [mm]')
# ax[0].set_xlabel('Time [sec]')
ax[0].legend(loc='lower right')
ax[0].grid()

ax1twin = ax[0].twinx()
ax1twin.plot(t, Cms_comp*1000, 'k', linewidth=1, label='Cms compensation')
ax1twin.axhline(y=Cms_min*1000, color='k', linestyle='--', linewidth=1, label='Lowest Cms compensation value')
#ax1twin.plot(t, track_Rms, 'k', linewidth=1, label='Rms compensation')
ax1twin.set_ylim([Cms_min*990, Cms*1100])
ax1twin.legend(loc='upper right')
ax1twin.set(ylabel='Compliance [mm/N]')
# ax1twin.set(ylabel='Rms comp [kg/s]')


# ax[1].plot(t, x2, linewidth=1, alpha=0.2, label='Displacement (no feedback)')
ax[1].plot(t, rms, 'm', linewidth=2, alpha=0.6, label=f'RMS estimator (no feedback)\n Tavg {np.int16(T_avg*1000)} ms')
ax[1].axhline(y=Threshold, color='m', linestyle='--', alpha=0.6, linewidth=1, label=f'Threshold: {Threshold} mm')
ax[1].set_ylabel('Amplitude [mm]')
ax[1].set_xlabel('Time [sec]')
ax[1].legend(loc='center right')
ax[1].grid()

ax1twin = ax[1].twinx()
ax1twin.plot(t, track_AutoR, label='Auto release ratio (feedback)')
ax1twin.set(ylabel='Auto release ratio [none]')
ax1twin.legend(loc='upper right')

# plt.savefig(f'Figures/physical_approach/simulation_{loudspeaker}_settings3.pdf')
plt.show()