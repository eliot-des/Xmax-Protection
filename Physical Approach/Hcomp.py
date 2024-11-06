import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import os
import sys
fpath = os.path.join(os.getcwd(), 'Modules')
sys.path.append(fpath)
from filters import bilinear2ndOrder

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


#================================================================================
#Test gain match

Xmax_ = 1e-3
A = 6

R = (Xmax_*Rec/(A*Bl*Cms))
dBgain = 20*np.log10(R)
print(round(dBgain,2))

if R <= 1:
    pass
else:
    raise ValueError('Gain is too low or Xmax criteria always satisfied.')


#================================================================================
#Compensation filter

f   = np.geomspace(10, 20000, 100)
w   = 2*np.pi*f
Fs = 48000

Cms_comp = R*Cms

Q0   = 1/np.sqrt(2)                                  # neutral quality factor
Q0_c = 1/(Rms+(Bl**2)/Rec)*np.sqrt(Mms/Cms_comp)     # Quality factor for the compensation filter - c stands for compensation

if Q0_c > Q0:
        Rms_comp = 1/Q0*np.sqrt(Mms/Cms_comp) - (Bl**2/Rec)
else:
        Rms_comp = Rms

b = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])         
a = np.array([Mms, Rms_comp+Bl**2/Rec, 1/Cms_comp])

bd, ad = bilinear2ndOrder(b, a, Fs)

_, H = sig.freqz(bd, ad, worN=f, fs=Fs)



#================================================================================
#Plot it

fig, ax = plt.subplots()
ax.semilogx(f, 20*np.log10(np.abs(H)), 'b', label='H_comp')
ax.axhline(y=dBgain, color='red', linestyle='--', label='Gain from Cms ratio')
ax.set_xlim([f[0], f[-1]])
ax.set_xlabel('Frequency  [Hz]')
ax.set_ylabel('Gain [dB]')
ax.grid(which='both')
ax.legend(loc='lower right')
plt.show()