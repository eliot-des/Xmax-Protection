import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
from lib import normalize, multibandChebyshev1Filters, multibandEllipticFilters, dynamic_peak_follower1, dynamic_rms_follower

#================================================================================
# The following code reads an audio file and processes it with a multi-band.
# The audio file is considered to be a voltage signal send to a loudspeaker.
# The displacement of the loudspeaker is calculated using Thiele-Small parameters
# by applying a digital filter to the discrete voltage signal.

# The displacement signal is then processed with a multi-band filter.
# The peak of the displacement signal is calculated for each band.
#================================================================================

Fs, v = read('Audio Tests/Thriller.wav')
v = v[:,0]
v = normalize(v)    #normalize the signal to 1

G = 10               #gain of the amplifier -> Max tension in volts
u = G*v             #tension in volts

t = np.arange(0, len(v)/Fs, 1/Fs)
tstart   = 0            #start time in seconds
duration = 1            #duration time in seconds

u = u[int(tstart*Fs):int((tstart+duration)*Fs)]
t = t[int(tstart*Fs):int((tstart+duration)*Fs)]
#================================================================================

#Defining Xmax
Xmax = 0.0012 

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

#Displacement filter
f   = np.geomspace(1, Fs/2, 1000)
w   = 2*np.pi*f
Hxu = Bl/((Rec+1j*w*Lec)*(-(w**2)*Mms+1j*w*Rms+1/Cms)+1j*w*Bl**2)

#analog coefficients
b = np.array([0, 0, 0, Bl])
a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])

#convert to digital filter
b, a = sig.bilinear(b, a, Fs)

x = sig.lfilter(b, a, u)    #displacement in m
x*=1e3                     #displacement in mm

#================================================================================
# Multi-band processing
fc = [40, 90, 200]     #cutoff frequencies

attack_time  = 0.00005      # Attack time in seconds
release_time = 0.07         # Release time in seconds
averaging_time = 0.05       # Averaging time in seconds

#sos_filters = multibandChebyshev1Filters(4, fc, Fs, 1)
sos_filters = multibandEllipticFilters(4, fc, Fs, 1, 60)

x_filt      = np.zeros((len(sos_filters), len(x)))
x_peak      = np.zeros((len(sos_filters), len(x)))
x_rms       = np.zeros((len(sos_filters), len(x)))
H_chebyshev = np.zeros((len(sos_filters), len(f)), dtype=complex)

for i in range(len(sos_filters)):
    x_filt[i] = sig.sosfilt(sos_filters[i], x)
    x_peak[i] = dynamic_peak_follower1(x_filt[i], attack_time, release_time, Fs)
    x_rms[i]  = dynamic_rms_follower(x_filt[i], averaging_time, Fs)
    _, H_chebyshev[i] = sig.sosfreqz(sos_filters[i], worN = f, fs = Fs)
#================================================================================
# Plot

colors = plt.cm.Set2(np.arange(len(sos_filters)))
#create a subplot mosaic


fig, ax = plt.subplot_mosaic("AB;CD;EF", figsize=(10, 6), layout='constrained')

#set general title
fig.suptitle(f'Displacement prediction, peak follower for {loudspeaker} driver, for max voltage = {G} V.')

ax['A'].plot(t, u,'k', label='Tension',alpha=0.3)
ax['A'].set(xlabel='Time [s]', ylabel='Amplitude [V]', xlim=(tstart, tstart+duration))
ax['A'].grid()
ax['A'].legend(loc='upper left')
ax['A'].set_yticks(np.linspace(-G*1.02, G*1.02, 7))
ax['A'].set(ylim=(-G*1.02, G*1.02))

Axtwin = ax['A'].twinx()
Axtwin.plot(t, x, label='Displacement')
Axtwin.plot(t, Xmax*np.ones(len(t))*1e3, 'r--',alpha=0.5, label='Xmax')
Axtwin.set(ylabel='Displacement [mm]')
Axtwin.legend(loc='lower left')

a = max(-Axtwin.get_yticks()[0], Axtwin.get_yticks()[-1])
Axtwin.set_yticks(np.linspace(-a, a, len(ax['A'].get_yticks())))



for i in range(len(sos_filters)):
    ax['B'].semilogx(f, 20*np.log10(np.abs(H_chebyshev[i])), label=f'Band {i+1}', color = colors[i])
    ax['B'].set(xlabel='Frequency [Hz]', ylabel='Magnitude [dB]')
    ax['B'].grid(True, which='both')
    ax['B'].legend(loc='upper right')
    ax['B'].set(xlim=(20, 20000), ylim=(-45, 5))


Bxtwin = ax['B'].twinx()
Bxtwin.semilogx(f, np.abs(Hxu)*1e3, 'k', alpha = 0.1)
Bxtwin.set(ylabel='|X/U| [mm/V]')

for i in range(len(sos_filters)):
    letter = chr(ord('C') + i)

    ax[letter].plot(t, x_filt[i], label='$x_{{filt}}^{}$'.format(i+1), color = colors[i])
    ax[letter].plot(t, x_peak[i], 'k', label='Peak($x_{{filt}}^{}$)'.format(i+1))
    ax[letter].plot(t, x_rms[i], color='crimson', label='RMS($x_{{filt}}^{}$)'.format(i+1))
    ax[letter].set(xlabel='Time [s]', ylabel='Displacement [mm]')
    ax[letter].grid()
    ax[letter].legend()
    ax[letter].set(xlim=(tstart, tstart+duration))
plt.show()

