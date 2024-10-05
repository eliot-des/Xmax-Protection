import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize, dynamic_peak_follower2, gain_factor_smoothing_sbs_bis, dynamic_peak_follower_sbs
import sys

plt.rc('lines', linewidth=2)
plt.rc('font', size=14)
plt.rc('axes', linewidth=1.5, labelsize=14)
plt.rc('legend', fontsize=14)

#================================================================================
#============= Definition of a signal being the input voltage ===================
#================================================================================
'''
Fs, v = read('Audio Tests/Thriller.wav')
t = np.arange(0, len(v)/Fs, 1/Fs)
'''
#create a sweep/ chirp signal
Fs = 48000
t = np.arange(0, 1, 1/Fs)
v = sig.chirp(t, f0=10000, f1=10, t1=1, method='logarithmic')
#v = normalize(v)    #normalize the signal to 1

A = 4               #gain of the amplifier -> Max tension in volts
u = A*v             #tension in volts

# fig, ax = plt.subplots()
# ax.plot(t, u, label='Sweep')
# ax.set_xlabel('Time [sec]')
# ax.set_ylabel('Amplitude [u.a]')
# ax.legend()
# plt.show()

# sys.exit()


tstart   = 0            #start time in seconds
duration = 1            #duration time in seconds

u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]


#================================================================================
#==== Importing speaker TS parameters and defining the displacement filter ======
#================================================================================

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

#Displacement calculation without feedback control

# x2 = sig.lfilter(b, a, u)    #displacement in m

x2 = np.zeros_like(u)
d = np.zeros(4) # memory buffer for Direct Form I

for n in range(len(u)):
    x2[n] = b[0]*u[n] + b[1]*d[0] + b[2]*d[1] - a[1]*d[2] - a[2]*d[3]

    d[1] = d[0]
    d[0] = u[n]
    d[3] = d[2]
    d[2] = x2[n]

x2*= 1e3                     #displacement in mm


# fig, ax = plt.subplots()
# ax.plot(t, u, label='Sweep')
# ax.plot(t, x2, label='Estimated displacement')
# ax.set_xlabel('Time [sec]')
# ax.set_ylabel('Amplitude [a.u]')
# ax.legend()
# ax.grid()
# plt.show()

# sys.exit()


#================================================================================
#====================== Compensation filter definition ==========================
#========================= Initializing parameters ==============================
#================================================================================

#feedback control algorithm

#Initializations

Xmax = 1 #in mm 
Xmargin = 0.5

x        = np.zeros(len(u))    #displacement in mm
x_lim    = np.zeros(len(u))    #limited displacement in mm
x_peak   = np.zeros(len(u))    #enveloppe estimator in mm
Cms_comp = np.zeros(len(u))    #compensation compliance to be adjusted sample by sample

Cms_max = 0.95*Cms             #initial compensation compliance equal (almost) to the real one of the speaker
Cms_min = 0.4*Cms              #minimum compensation compliance is 10% of the real Cms (CHOICE)
Cms_comp[0] = Cms_max

attack_peak    = 0.00005     # Attack time in seconds
release_peak   = 0.1       # Release time in seconds
attack_smooth  = 0.002     # Attack time for the gain smoothing function
release_smooth = 0.005      # Release time for the gain smoothing function

# print("k attack: {} sec".format(round(1 - np.exp(-2.2 / (attack_smooth * Fs)), 3)))
# print("k release: {} sec".format(round(1 - np.exp(-2.2 / (release_smooth * Fs)), 3)))

# sys.exit()

#create the compensation filter based on the physical approach
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
# ax.legend()
# plt.show()

# sys.exit()


#================================================================================
#================================== Main loop ===================================
#================================================================================

d_DFI = np.zeros(4)     # memory buffer for Direct Form I
d_TDFII = np.zeros(2)   # memory buffer for Transposed Direct Form I
a_comp = np.zeros(3)    # denumerator coefficient of the compensation filter in the z domain

for i in range(1, len(u)):
    
    #we first estimate the displacement from the u[i] sample
    x[i] = b[0]*u[i] + b[1]*d_DFI[0] + b[2]*d_DFI[1] - a[1]*d_DFI[2] - a[2]*d_DFI[3]

    d_DFI[1] = d_DFI[0]
    d_DFI[0] = u[i]
    d_DFI[3] = d_DFI[2]
    d_DFI[2] = x[i]

    x_peak[i] = dynamic_peak_follower_sbs(x[i]*1000, x_peak[i-1], attack_peak, release_peak, Fs)
    
    #we then compare the peak of the displacement signal with Xmax
    # if np.abs(x[i]*1000) > Xmax:
    if np.abs(x_peak[i]) > Xmax:
        # Cms_comp[i] = Cms_min
        Cms_target = Cms_min
    else:
        # Cms_comp[i] = Cms_max
        Cms_target = Cms_max

    # we then apply the gain smoothing function to the Cms_comp value of the compensation filter
    Cms_comp[i] = gain_factor_smoothing_sbs_bis(Cms_target, Cms_comp[i-1], attack_smooth, release_smooth, Fs)

    # #we then compute the new compensation coeff filter (zeros are unchanged)...
    A_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms_comp[i]])

    #convert to digital filter
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

    # #... then apply it to the tension signal
    x_lim[i] = b_comp[0]*x[i] + d_TDFII[0]

    d_TDFII[0] = b_comp[1]*x[i] - a_comp[1]*x_lim[i] + d_TDFII[1]
    d_TDFII[1] = b_comp[2]*x[i] - a_comp[2]*x_lim[i]

# sys.exit()


#================================================================================
#================================== Main plot ===================================
#================================================================================

fig, ax = plt.subplots()

ax.plot(t, x*1000, '-.b', linewidth=1, label='Estimated displacement (feedback)')
ax.plot(t, x_peak, 'r-', linewidth=1, alpha=0.4, label='x_{peak}')
ax.plot(t, x_lim*1000, 'm', label='Limited displacement (feedback)')
ax.axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
ax.axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
ax.set_ylabel('Amplitude [mm]')
ax.set_xlabel('Time [sec]')
ax.legend(loc='lower left')
ax.grid()

ax1twin = ax.twinx()
ax1twin.plot(t, Cms_comp, 'k', linewidth=1, label='Cms compensation')
ax1twin.legend(loc='upper left')
ax1twin.set(ylabel='Compliance [m/N]')

plt.show()

sys.exit()

#================================================================================
#================================ Main plot bis =================================
#================================================================================

fig, ax = plt.subplots(2, 1, sharex = True)

ax[0].plot(t, u, label='Tension | Feed.')
ax[0].set(xlabel='Time [s]', ylabel='Amplitude [V]', xlim=(tstart, tstart+duration))

ax[0].set_yticks(np.linspace(-A*1.02, A*1.02, 7))
ax[0].grid()

ax[1].plot(t, x2,'k', label='Displacement | No Feed.')
ax[1].plot(t, x,label='Displacement | Feed.')
# ax[1].plot(t, x_peak, 'r--',alpha=0.5, label='x_{peak}')
ax[1].set(xlabel='Time [s]', ylabel='Displacement [mm]', xlim=(tstart, tstart+duration))
ax[1].grid()
ax[1].legend()

ax1twin = ax[1].twinx()
ax1twin.plot(t, Cms_comp,'r', label='Cms compensation')
ax1twin.set(ylabel='Compliance [m/N]')

plt.show()