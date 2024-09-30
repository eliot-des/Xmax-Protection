import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize, dynamic_peak_follower_sbs, gain_factor_smoothing_sbs, PeakFilterTD2, HighPassFilterTD2

#================================================================================
'''
Fs, v = read('Audio Tests/Thriller.wav')
t = np.arange(0, len(v)/Fs, 1/Fs)
'''
#create a sweep/ chirp signal
Fs = 48000
t = np.arange(0, 1, 1/Fs)
v = sig.chirp(t, f0=20000, f1=10, t1=1, method='logarithmic')
#v = normalize(v)    #normalize the signal to 1

A = 10               #gain of the amplifier -> Max tension in volts
u = A*v             #tension in volts


tstart   = 0            #start time in seconds
duration = 1            #duration time in seconds

u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]
#================================================================================

#Defining Xmax and margin
Xmax = 2 #in mm 
Xmargin = 1

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
#Displacement calculation without feedback control

x2 = sig.lfilter(b, a, u)    #displacement in m
x2*= 1e3                     #displacement in mm

#================================================================================

#feedback control algorithm

#Initializations
x       = np.zeros(len(u))    #displacement in mm
x_peak  = np.zeros(len(u))    #peak of the displacement in mm
fc      = np.zeros(len(u))    #cut-off frequency in Hz of the low-pass filter

fc[0] = 10                    #initial cut-off frequency in Hz
fc_max = 200                  #maximum cut-off frequency in Hz


attack_peak    = 0.00005    # Attack time in seconds
release_peak   = 0.07       # Release time in seconds
attack_smooth  = 0.0001     # Attack time for the gain smoothing function
release_smooth = 0.005      # Release time for the gain smoothing function

#biquad function that will be use for the peak filter
#highPassFilter = HighPassFilterTD2(Fc=fc[0],  Q=0.72, dBgain=0, fs=Fs)
#create an high-pass filter with the scipy
b_hp, a_hp = sig.butter(2, 10, 'hp', fs=Fs)

f = np.geomspace(1, Fs/2, 1000)

fig, ax = plt.subplots()

_, h = sig.freqz(b_hp, a_hp, worN=f, fs=Fs)
ax.semilogx(f, 20*np.log10(np.abs(h)))
ax.grid(which='both')
plt.show()



ax.plot()


zi = sig.lfilter_zi(b, a) * u[0]            # Scale by the first sample of u to initialize properly
zi_hp = sig.lfilter_zi(b_hp, a_hp) * u[0]   # Scale by the first sample of u to initialize properly

#Main loop
for i in range(1, len(u)):
    
    #we first estimate the displacement from the u[i] sample
    x[i], zi = sig.lfilter(b, a, [u[i-1]], zi=zi)    #displacement in m
    x[i]*=1e3                                        #displacement in mm

    #x_peak[i] = dynamic_peak_follower_sbs(x[i], x_peak[i-1], attack_peak, release_peak, Fs)
    
    #we then compare the peak of the displacement signal with Xmargin
    if np.abs(x[i]) > Xmargin:
        fc_target = fc_max
    else:
        fc_target = 10

    #we then apply the gain smoothing function to the cut-off frequency of the high-pass filter
    fc[i] = gain_factor_smoothing_sbs(fc_target, fc[i-1], attack_smooth, release_smooth, Fs)

    #we then apply the high-pass filter with the new cut-off frequency to the tension signal
    u[i], zi_hp = sig.lfilter(b_hp, a_hp, [u[i]], zi=zi_hp)



fig, ax = plt.subplots(2, 1, sharex = True)

ax[0].plot(t, u, label='Tension | Feed.')
ax[0].set(xlabel='Time [s]', ylabel='Amplitude [V]', xlim=(tstart, tstart+duration))

ax[0].set_yticks(np.linspace(-A*1.02, A*1.02, 7))
ax[0].grid()

ax[1].plot(t, x2,'k', label='Displacement | No Feed.')
ax[1].plot(t, x,label='Displacement | Feed.')
ax[1].plot(t, x_peak, 'r--',alpha=0.5, label='x_{peak}')
ax[1].set(xlabel='Time [s]', ylabel='Displacement [mm]', xlim=(tstart, tstart+duration))
ax[1].grid()
ax[1].legend()

ax1twin = ax[1].twinx()
ax1twin.plot(t, fc,'r', label='Cut-off frequency')
ax1twin.set(ylabel='Frequency [Hz]')

plt.show()

