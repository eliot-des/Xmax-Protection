import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize, dynamic_peak_follower2, gain_factor_smoothing_sbs_bis
import sys

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

#Defining Xmax and margin
Xmax = 1 #in mm 
Xmargin = 0.5

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
b = np.array([0, 0, 0, Bl])
a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])

#convert to digital filter
b, a = sig.bilinear(b, a, Fs)


#================================================================================
#========== Filter the input voltage to get an estimated displacment ============
#================================================================================

#Displacement calculation without feedback control

x2 = sig.lfilter(b, a, u)    #displacement in m
x2*= 1e3                     #displacement in mm

# fig, ax = plt.subplots()
# ax.plot(t, u, label='Sweep')
# ax.plot(t, x2, label='Estimated displacement')
# ax.axhline(y=Xmax, color='red', linestyle='--', label='Xmax')
# ax.set_xlabel('Time [sec]')
# ax.set_ylabel('Amplitude [a.u]')
# ax.legend()
# plt.show()

# sys.exit()


#================================================================================
#====================== Compensation filter definition ==========================
#========================= Initializing parameters ==============================
#================================================================================

#feedback control algorithm

#Initializations
x        = np.zeros(len(u))    #displacement in mm
x_lim   = np.zeros(len(u))     #limited displacement in mm
Cms_comp = np.zeros(len(u))    #compensation compliance to be adjusted sample by sample

Cms_max = 0.95*Cms        #initial compensation compliance equal (almost) to the real one of the speaker
Cms_min = 0.1*Cms             #minimum compensation compliance is 10% of the real Cms (CHOICE)
Cms_comp[0] = Cms_max

# attack_peak    = 0.00005     # Attack time in seconds
# release_peak   = 0.07       # Release time in seconds
attack_smooth  = 0.0005     # Attack time for the gain smoothing function
release_smooth = 0.0005      # Release time for the gain smoothing function

# print("k attack: {} sec".format(round(1 - np.exp(-2.2 / (attack_smooth * Fs)), 3)))
# print("k release: {} sec".format(round(1 - np.exp(-2.2 / (release_smooth * Fs)), 3)))

# sys.exit()

#create the compensation filter based on the physical approach
b_comp = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])
a_comp = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms_comp[0]+Bl**2, Rec/Cms_comp[0]])

f = np.geomspace(1, Fs/2, 1000)
w   = 2*np.pi*f

# _, h = sig.freqs(b_comp, a_comp, worN=w)

# fig, ax = plt.subplots()

# ax.semilogx(f, 20*np.log10(np.abs(h)))
# ax.set_xlabel('Frequency  [Hz]')
# ax.set_ylabel('Gain [dB]')
# ax.grid(which='both')
# plt.show()

# sys.exit()


#================================================================================
#================================== Main loop ===================================
#================================================================================

zi = sig.lfilter_zi(b, a) * u[0]            # Scale by the first sample of u to initialize properly
zi_comp = sig.lfilter_zi(b_comp, a_comp) * u[0]   # Scale by the first sample of u to initialize properly

# sys.exit()

for i in range(1, len(u)):
    
    #we first estimate the displacement from the u[i] sample
    x[i], zi = sig.lfilter(b, a, [u[i-1]], zi=zi)    #displacement in m
    x[i]*=1e3                                        #displacement in mm

    #x_peak[i] = dynamic_peak_follower_sbs(x[i], x_peak[i-1], attack_peak, release_peak, Fs)
    
    #we then compare the peak of the displacement signal with Xmax
    if np.abs(x[i]) > Xmax:
        # Cms_comp[i] = Cms_min
        Cms_target = Cms_min
    else:
        # Cms_comp[i] = Cms_max
        Cms_target = Cms_max

    # we then apply the gain smoothing function to the Cms_comp value of the compensation filter
    Cms_comp[i] = gain_factor_smoothing_sbs_bis(Cms_target, Cms_comp[i-1], attack_smooth, release_smooth, Fs)

    # #we then compute the new compensation coeff filter (zeros are unchanged)...
    a_comp = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms_comp[i]+Bl**2, Rec/Cms_comp[i]])

    #convert to digital filter
    # b_comp_d, a_comp_d = sig.bilinear(b_comp, a_comp, Fs)

    # #... then apply it to the tension signal
    # x_lim[i], zi_comp = sig.lfilter(b_comp_d, a_comp_d, [u[i]], zi=zi_comp)


# sys.exit()

fig, ax = plt.subplots(2, 1)

ax[0].plot(t, u, label='Sweep')
# ax[0].plot(t, x2+0.05, label='Estimated displacement (whole array)')
ax[0].plot(t, x, label='Estimated displacement (feedback)')
ax[0].axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
ax[0].axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
ax[0].set_ylabel('Amplitude [a.u]')
ax[0].legend(loc='lower left')

ax[1].plot(t, x_lim*1000, label='Limited displacement')
ax[1].set_ylim([-10, 10])
ax[1].set_xlabel('Time [sec]')
ax[0].set_ylabel('Amplitude [mm]')

ax1twin = ax[0].twinx()
ax1twin.plot(t, Cms_comp, 'r', label='Cms compensation')
ax1twin.legend(loc='upper left')
ax1twin.set(ylabel='Compliance [m/N]')

plt.show()

sys.exit()

#================================================================================
#================================== Main plot ===================================
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
# %%
