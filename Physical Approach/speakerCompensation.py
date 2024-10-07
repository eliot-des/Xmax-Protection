import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

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

#Displacement filter
f   = np.geomspace(1, 20000, 1000)
w   = 2*np.pi*f
Hxu = Bl/((Rec+1j*w*Lec)*(-(w**2)*Mms+1j*w*Rms+1/Cms)+1j*w*Bl**2)

#analog coefficients
b = np.array([0, 0, 0, Bl])
a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])



#compensation analog filter
Xmax = 1.2e-3

w0 = 1/np.sqrt(Mms*Cms)
Q0 = 1/(w0*Cms*(Rms+(Bl**2/Rec)))
#Q0 = 1/np.sqrt(2)
A_list = np.geomspace(10, 30, 10)


Mms_comp = Mms

b_comp = a
a_comp = [None]*len(A_list)

""
for i in range(len(A_list)):
    
    Cms_comp = Xmax*Rec/(A_list[i]*Bl)
    Rms_comp = 1/Q0 * np.sqrt(Mms/Cms_comp) - (Bl**2/Rec)

    a_comp[i] = [Lec*Mms_comp, Rec*Mms_comp+Lec*Rms_comp, Rec*Rms_comp+Lec/Cms_comp+Bl**2, Rec/Cms_comp]


#frequency response

_, h = sig.freqs(b, a, worN=w)


h_comp = [None]*len(A_list)

for i in range(len(A_list)):
    _, h_comp[i] = sig.freqs(b_comp, a_comp[i], worN=w)



fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].semilogx(f, np.abs(h)*1000, '--r', label='Original')

#create color map for the bands
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(A_list)))

for i in range(len(A_list)):
    ax[0].semilogx(f, np.abs(h*h_comp[i])*1000, color = colors[i], label=f'Cms: {A_list[i]:.2e} m/N')
    ax[1].semilogx(f, np.rad2deg(np.angle(h_comp[i])), color = colors[i])

ax[0].set(xlim = (20, 210)) #, ylim = (-45, 3)
ax[0].set_ylabel('Magnitude [mm/V]')
# ax[0].set_title()
ax[0].grid(which='both', axis='both')
ax[0].legend(loc='upper right')

ax[1].set(xlabel = 'Frequency [Hz]', ylabel = 'Angle [deg]')
ax[1].set(xlim = (10, f[-1])) #, ylim = (0, 100)
ax[1].grid(which='both', axis='both')

plt.savefig(f'Figures/compensation/compensation_{loudspeaker}.pdf')
plt.show()