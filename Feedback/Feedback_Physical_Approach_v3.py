import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read, write
import os
import sys
fpath = os.path.join(os.getcwd(), 'modules')
sys.path.append(fpath)
from dynamic_trackers import dynamic_peak_follower2, gain_factor_smoothing_sbs_bis, dynamic_peak_follower_sbs, gain_factor_smoothing_sbs
from audio import normalize
from filters import bilinear2ndOrder

plt.rc('lines', linewidth=2)
plt.rc('font', size=14)
plt.rc('axes', linewidth=1.5, labelsize=14)
plt.rc('legend', fontsize=12)

'''
Feddback V3 dot not consider resonant speaker. Compensation filter quality
factor is adjusted to stay under 1/sqrt(2)
'''


#================================================================================
#============= Definition of a signal being the input voltage ===================
#================================================================================

# signal_name = 'ShookOnes'
# Fs, u = read('Audio/'+signal_name+'.wav')
# t = np.arange(0, len(u[:,0])/Fs, 1/Fs)
# u = normalize(u[:,0])

# create a sweep/chirp signal
signal_name = 'S'   # to further name the .wav export file - S stands for sweep
Fs = 48000
t = np.arange(0, 5, 1/Fs)
u = sig.chirp(t, f0=10000, f1=20, t1=5, method='logarithmic')

A = 10                   # gain of the amplifier -> Max tension in volts
u = A*u                 # tension in volts
 
# tstart   = 11            # start time in seconds
# duration = 3            # duration time in seconds

# u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
# t = t[int(tstart*Fs):int((tstart+duration)*Fs)]

#================================================================================
#==== Importing speaker TS parameters and defining the displacement filter ======
#================================================================================

Xmax = 1            #in mm
Xmax_= Xmax*1e-3    #in m

# Thiele-small parameters
loudspeakers = {'full-range1' :'Dayton_CE4895-8',
                'full-range2' :'Dayton_HARB252-8',
                'woofer1': 'Peerless_HDSP830860',
                'woofer1_Klippel': 'Peerless_HDSP830860_Klippel',
                'woofer2': 'Dayton_DCS165-4',
                'woofer3': 'Dayton_RS150-4',
                'subwoofer1': 'B&C_15FW76-4'}

loudspeaker = loudspeakers['woofer1_Klippel']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)

# print(Rms)
# sys.exit()

# Low-frequency approximation of displacement
B_LF = np.array([0, 0, Bl/Rec])         
A_LF = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])


# Analytical z-domain coeff from Low-frequency approximation

b, a = bilinear2ndOrder(B_LF, A_LF, Fs)

#================================================================================
#========== Filter the input voltage to get an estimated displacement ===========
#================================================================================
x2 = np.zeros_like(u)        # displacement in m
x2 = sig.lfilter(b, a, u)    # displacement in m
x2 = x2*1000                 # displacement in mm

#================================================================================
#====================== Compensation filter definition ==========================
#========================= Initializing parameters ==============================
#================================================================================

C = Xmax_*Rec/(A*Bl*Cms)            # optimal minimal compliance to respect Xmax
print('Chosen C = ', np.round(C, 3))

if C <= 1:
    pass
else:
    raise ValueError('Gain is too low or Xmax criteria always satisfied.')

Cms_comp = np.zeros(len(u))         # compensation compliance to be adjusted sample by sample
Cms_comp[0] = Cms                   # initial compliance value is the one of the speaker
Cms_min = 0.9*C*Cms                 # minimum compensation compliance with a margin factor

Rms_comp = Rms*np.ones_like(u)      # compensation loss to be adjusted sample by sample


B_comp = A_LF
A_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms_comp[0]])

b_comp, a_comp = bilinear2ndOrder(B_comp, A_comp, Fs)

Q0 = 1/np.sqrt(2)                               # neutral quality factor
Q0_s = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms)     # original loudspeaker quality factor - s stands for speaker

# print("Speaker quality factor:", round(Q0_s,3))
# print("Comp filter numerator coeff:"  , b_comp)
# print("Comp filter denumerator coeff:", a_comp)
# sys.exit()


#================================================================================
#================== Main loop - feedback control algorithm ======================
#================================================================================

#Initializations

x        = np.zeros(len(u))    # displacement in mm

attack_smooth  = 0.02          # Attack time for the gain smoothing function
release_smooth = 0.5           # Release time for the gain smoothing function

d_DFI   = np.zeros(4)          # memory buffer for Direct Form I
d_TDFII = np.zeros(2)          # memory buffer for Transposed Direct Form I
u_hp    = np.zeros_like(u)

for i in range(1, len(u_hp)):
    
    #we first estimate the displacement from the u_hp[i] sample
    x[i] = b[0]*u_hp[i-1] + b[1]*d_DFI[0] + b[2]*d_DFI[1] - a[1]*d_DFI[2] - a[2]*d_DFI[3]

    d_DFI[1] = d_DFI[0]
    d_DFI[0] = u_hp[i-1]
    d_DFI[3] = d_DFI[2]
    d_DFI[2] = x[i]

    x[i] *= 1000

    #we then compare the peak of the displacement signal with Xmax
    if np.abs(x[i]) > Xmax:
        Cms_target = Cms_min
    else:
        Cms_target = Cms

    # we then apply the gain smoothing function to the Cms_comp value of the compensation filter
    Cms_comp[i] = gain_factor_smoothing_sbs_bis(Cms_target, Cms_comp[i-1], attack_smooth, release_smooth, Fs)

    Q0_c = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms_comp[i])     # Quality factor for the compensation filter - c stands for compensation

    if Q0_c > Q0:
        Rms_comp[i] = 1/Q0*np.sqrt(Mms/Cms_comp[i]) - (Bl**2/Rec)
    else:
        Rms_comp[i] = Rms
    

    #we then compute the new compensation coeff filter (zeros are unchanged)...
    A_comp = np.array([Mms, Rms_comp[i]+Bl**2/Rec, 1/Cms_comp[i]])

    #compute the digital coefficients
    b_comp, a_comp = bilinear2ndOrder(B_comp, A_comp, Fs)
    

    #Apply the compensation filter to the input tension signal
    u_hp[i] = b_comp[0]*u[i] + d_TDFII[0]

    d_TDFII[0] = b_comp[1]*u[i] - a_comp[1]*u_hp[i] + d_TDFII[1]
    d_TDFII[1] = b_comp[2]*u[i] - a_comp[2]*u_hp[i]


#================================================================================
#=================================== Main plot ==================================
#================================================================================

fig, ax = plt.subplots(4, 1, sharex=True)

ax[0].plot(t, u, 'b', alpha=0.2, linewidth=1, label='U in')
ax[0].set_ylabel('[mm]')
ax[0].legend(loc='lower left')
ax[0].grid()

ax[1].plot(t, u_hp, linewidth=1, label='U hp')
ax[1].set_ylabel('[mm]')
ax[1].legend(loc='lower left')
ax[1].grid()

ax[2].plot(t, x, 'b', linewidth=1, label=f'Estimated displacement (feedback)\n Attack {attack_smooth*1e3:.0f} ms\n Release {release_smooth*1e3:.0f} ms')
ax[2].plot(t, x2, linewidth=1, alpha=0.2, label='Displacement (no feedback)')
ax[2].axhline(y=Xmax, color='red', linestyle='--', label='+/- Xmax')
ax[2].axhline(y=-Xmax, color='red', linestyle='--')
ax[2].set_ylabel('[mm]')
ax[2].legend(loc='lower left')
ax[2].grid()

ax[3].plot(t, Cms_comp/Cms, 'r', linewidth=1, label='Cms_comp/Cms')
ax[3].set_ylabel('C ratio')
ax[3].set_ylim([-0.1, 1.1])
ax[3].legend(loc='lower left')
ax[3].grid()

ax1twin = ax[3].twinx()
ax1twin.plot(t, Rms/Rms_comp, 'k', linewidth=1, label='Rms/Rms_comp')
ax1twin.set_ylim([-0.1, 1.1])
ax1twin.legend(loc='upper left')
ax1twin.set(ylabel='R ratio')

plt.show()


#================================================================================
#=========================== Main plot - old version ============================
#================================================================================

# fig, ax = plt.subplots(2,1, sharex=True)

# ax[0].plot(t, x, 'b', linewidth=1, label=f'Estimated displacement (feedback)\n Attack {attack_smooth*1e3:.0f} ms\n Release {release_smooth*1e3:.0f} ms')
# ax[0].plot(t, x2, linewidth=1, alpha=0.2, label='Displacement (no feedback)')
# ax[0].axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
# ax[0].axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
# ax[0].set_ylabel('Amplitude [mm]')
# # ax[0].legend(loc='lower left')
# ax[0].grid()

# ax1twin = ax[0].twinx()
# ax1twin.plot(t, Cms_comp*1000, 'k', linewidth=1, label='Cms compensation')
# # ax1twin.legend(loc='upper left')
# ax1twin.set(ylabel='Compliance [mm/N]')


# ax[1].plot(t, x, 'b', linewidth=1, label=f'Estimated displacement (feedback)\n Attack {attack_smooth*1e3:.0f} ms\n Release {release_smooth*1e3:.0f} ms')
# ax[1].plot(t, x2, linewidth=1, alpha=0.2, label='Displacement (no feedback)')
# ax[1].axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
# ax[1].axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
# ax[1].set_ylabel('Amplitude [mm]')
# ax[1].set_xlabel('Time [sec]')
# # ax[1].legend(loc='lower left')
# ax[1].grid()

# ax1twin = ax[1].twinx()
# ax1twin.plot(t, Rms_comp, 'k', linewidth=1, label='Rms compensation')
# ax1twin.axhline(y=Rms, color='m', linestyle='-.', label='Rms speaker')
# # ax1twin.legend(loc='upper left')
# ax1twin.set(ylabel='Rms comp [kg/s]')

# # plt.savefig(f'Figures/physical_approach/simulation_{loudspeaker}.pdf')
# plt.show()


#================================================================================
#============================== Writting data ===================================
#================================================================================

# Max_Uin = np.max(np.abs(u))
# Max_Uhp = np.max(np.abs(u_hp))

# norm_factor = np.maximum(Max_Uin, Max_Uhp)
# print('Normalization factor:', np.round(norm_factor, 1))
# print('Max u_in:', round(Max_Uin, 1))
# print('Max u_hp:', round(Max_Uhp, 1))

# u=u/(norm_factor)
# u_hp=u_hp/(norm_factor)

# u*=2**15
# u_hp*=2**15

# write(f'Audio/Feedback/{signal_name}/G{A}_mono.wav', Fs, u.astype(np.int16))
# write(f'Audio/Feedback/{signal_name}/G{A}_Att{np.int16(attack_smooth*1000)}ms_{loudspeaker}.wav', Fs, u_hp.astype(np.int16))