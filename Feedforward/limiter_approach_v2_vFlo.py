import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig 
from scipy.io.wavfile import read
import os
import sys
fpath = os.path.join(os.getcwd(), 'modules')
sys.path.append(fpath)
from dynamic_trackers import dynamic_peak_follower_sbs, gain_factor_smoothing_sbs_bis
from audio import normalize
from filters import bilinear2ndOrder

plt.rc('lines', linewidth=2)
plt.rc('font', size=14)
plt.rc('axes', linewidth=1.5, labelsize=14)
plt.rc('legend', fontsize=12)


#================================================================================
#============= Definition of a signal being the input voltage ===================
#================================================================================

# create a sweep/chirp signal
Fs = 48000
t = np.arange(0, 5, 1/Fs)
u = sig.chirp(t, f0=10000, f1=10, t1=5, method='logarithmic')

A = 8                   # gain of the amplifier -> Max tension in volts
u = A*u                 # tension in volts


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

loudspeaker = loudspeakers['full-range2']

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
#====================== Compensation filter definition ==========================
#========================= Initializing parameters ==============================
#================================================================================
x        = np.zeros(len(u))       # displacement in mm
x_lim    = np.zeros(len(u))       # limited displacement in mm
Cms_comp = Cms*np.ones(len(u))    # compensation compliance to be adjusted sample by sample
d_DFI    = np.zeros(4)            # memory buffer for Direct Form I
d_TDFII  = np.zeros(2)            # memory buffer for Transposed Direct Form I
Q_track  = np.zeros(len(u))       # to track the quality factor
Q_track[0] = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms_comp[0])

B_comp = A_LF
A_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms_comp[0]])

b_comp, a_comp = bilinear2ndOrder(B_comp, A_comp, Fs)


#================================================================================
#==================== parallel branch for main analysis =========================
#================================================================================
for i in range(1, len(u)):

    #we first estimate the displacement from the u[i] sample
    x[i] = b[0]*u[i-1] + b[1]*d_DFI[0] + b[2]*d_DFI[1] - a[1]*d_DFI[2] - a[2]*d_DFI[3]

    d_DFI[1] = d_DFI[0]
    d_DFI[0] = u[i-1]
    d_DFI[3] = d_DFI[2]
    d_DFI[2] = x[i]

    x[i] *= 1000

    #compute the required Cms_comp
    if np.abs(x[i]) < Xmax:
        Cms_comp[i] = Cms
    else:
        Cms_comp[i] = Xmax/np.abs(x[i])*Cms

    # check the Q factor
    Q_track[i] = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms_comp[i])

    #we then compute the new compensation coeff filter (zeros are unchanged)...
    A_comp = np.array([Mms, Rms+Bl**2/Rec, 1/Cms_comp[i]])

    #compute the digital coefficients
    b_comp, a_comp = bilinear2ndOrder(B_comp, A_comp, Fs)
    

    #Apply the compensation filter to the input tension signal
    x_lim[i] = b_comp[0]*x[i] + d_TDFII[0]

    d_TDFII[0] = b_comp[1]*x[i] - a_comp[1]*x_lim[i] + d_TDFII[1]
    d_TDFII[1] = b_comp[2]*x[i] - a_comp[2]*x_lim[i]


#================================================================================
#================================== Main plot ===================================
#================================================================================
fig, ax = plt.subplots()

ax.plot(t, x, 'b', linewidth=1, alpha=0.4, label='Input displacement')
ax.plot(t, x_lim, 'b', linewidth=1, label='Output displacement')
ax.axhline(y=Xmax, color='red', linestyle='--', label='+Xmax')
ax.axhline(y=-Xmax, color='red', linestyle='--', label='-Xmax')
ax.set_ylabel('Amplitude [mm]')
ax.legend(loc='lower left')
ax.grid()

ax1twin = ax.twinx()
ax1twin.plot(t, Cms_comp*1000, 'k', linewidth=1, alpha=0.4, label='Cms compensation')
ax1twin.axhline(y=Cms*1000, color='k', alpha=0.5, linestyle='-.', label='Cms speaker')
ax1twin.legend(loc='upper left')
ax1twin.set(ylabel='Compliance [mm/N]')

# ax1twin = ax.twinx()
# ax1twin.plot(t, Q_track, 'm', linewidth=1, label='Quality factor')
# ax1twin.axhline(y=1/np.sqrt(2), color='m', alpha=0.5, linestyle='-.', label='1/sqrt(2)')
# ax1twin.legend(loc='upper left')

plt.show()