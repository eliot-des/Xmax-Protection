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

# fig, ax = plt.subplots()
# ax.plot(t, u, label='Signal')
# ax.set_xlabel('Time [sec]')
# ax.set_ylabel('Amplitude [u.a]')
# ax.legend()
# plt.show()

# sys.exit()


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

loudspeaker = loudspeakers['full-range1']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)


#optimal minimal compliance to not 
R_v1 = (Xmax_*Rec/(A*Bl*Cms))
print('Chosen R = ', np.round(R_v1, 3))

f_R = 20
omega_R = 2*np.pi*f_R

print((Bl*A/(Rec*Xmax_))**2-omega_R**2*(Rms**2+2*Rms*Bl**2/Rec+Bl**4/Rec**2)>0)

Kms1 = omega_R**2*Mms - np.sqrt((Bl*A/(Rec*Xmax_))**2-omega_R**2*(Rms**2+2*Rms*Bl**2/Rec+Bl**4/Rec**2))
Kms2 = omega_R**2*Mms + np.sqrt((Bl*A/(Rec*Xmax_))**2-omega_R**2*(Rms**2+2*Rms*Bl**2/Rec+Bl**4/Rec**2))
R = 1/(Cms*Kms2)

print(f'Root n°1 = {1/Kms1:.2e}, Root n°2 = {1/Kms2:.2e}')
print(f'Ratio n°1 = {np.round(1/(Cms*Kms2), 3)}, Ratio n°2 = {np.round(1/(Cms*Kms1), 3)}')
print('Chosen R = ', np.round(R, 3))

# sys.exit()


#analog coefficients
# b = np.array([0, 0, 0, Bl])
# a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])

B_LF = np.array([0, 0, Bl/Rec])         # Low-frequency approximation
A_LF = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])

#convert to digital filter
# b, a = sig.bilinear(b, a, Fs)

b = np.zeros(3)                         # Analytical z-domain coeff from Low-frequency approximation
b[0] = B_LF[2]/Fs**2
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
#========== Filter the input voltage to get an estimated displacment ============
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

# _, h = sig.freqs(B_comp, A_comp, worN=w)
# _, H = sig.freqz(b_comp, a_comp, worN=f, fs=Fs)

# fig, ax = plt.subplots()
# ax.semilogx(f, 20*np.log10(np.abs(h)), '--r', label='theoritical')
# ax.semilogx(f, 20*np.log10(np.abs(H)), '-.b', label='analytically bilinearized')
# ax.set_xlabel('Frequency  [Hz]')
# ax.set_ylabel('Gain [dB]')
# ax.grid(which='both')
# ax.legend(loc='lower right')
# plt.show()
# sys.exit()

x_lim = np.zeros_like(x2)
d_TDFII = np.zeros(2)   # memory buffer for Transposed Direct Form II

for i in range(len(x2)):
    x_lim[i] = b_comp[0]*x2[i] + d_TDFII[0]

    d_TDFII[0] = b_comp[1]*x2[i] - a_comp[1]*x_lim[i] + d_TDFII[1]
    d_TDFII[1] = b_comp[2]*x2[i] - a_comp[2]*x_lim[i]

# fig, ax = plt.subplots()
# ax.plot(t, u, label='Sweep')
# ax.plot(t, x2, label='Displacement')
# ax.plot(t, x_lim, label=f'Filtered displacement - Cms bis ={np.int16(ratio*100)}% Cms0')
# ax.set_xlabel('Time [sec]')
# ax.set_ylabel('Amplitude [a.u]')
# ax.legend(loc='lower left')
# ax.grid()
# plt.show()
# sys.exit()


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
# print("k attack: {} sec".format(round(1 - np.exp(-2.2 / (attack_smooth * Fs)), 3)))
# print("k release: {} sec".format(round(1 - np.exp(-2.2 / (release_smooth * Fs)), 3)))

# sys.exit()

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

    # x_peak[i] = dynamic_peak_follower_sbs(x[i], x_peak[i-1], attack_peak, release_peak, Fs) # x_peak in mm
    
    # Klippel trick
    dx = (x[i] - x[i-1]) * Fs
    x_peak[i] = np.sqrt(x[i]**2 + (dx/wr)**2)

    #we then compare the peak of the displacement signal with Xmax
    # if np.abs(x[i]) > Xmax:
    if x_peak[i] > Xmax:
        # Cms_comp[i] = Cms_min
        Cms_target = Cms_min
    else:
        # Cms_comp[i] = Cms_max
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

# sys.exit()


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