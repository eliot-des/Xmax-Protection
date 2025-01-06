import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import sys

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "legend.fontsize":15})

plt.close("all")


# =================== KLIPPEL MEASUREMENT==================
Z_klippel   = 'Measurements/Klippel/Magnitude_of_electric_impedance_Z(f).txt'
Hxu_Klippel = 'Measurements/Klippel/Magnitude_of_transfer_function_Hx(f).txt'

rawdata_Z = np.loadtxt(Z_klippel, skiprows=3)
rawdata_H = np.loadtxt(Hxu_Klippel, skiprows=5)

f_axis_1  = rawdata_Z[:,0]    # 1st column = Frequency (Hz)
Z_m       = rawdata_Z[:,1]    # 2nd column = Z measured (Ohm)
Z_f       = rawdata_Z[:,3]    # 4th column =  Z fitted (Ohm)

f_axis_2  = rawdata_H[:,0]    # 1st column = Frequency (Hz)
H_m       = rawdata_H[:,1]    # 2nd column = Measured (mm/V)
H_f_noC   = rawdata_H[:,3]    # 4th column =  Z fitted without crepp (mm/V)
H_f       = rawdata_H[:,5]    # 5th column =  Z fitted (mm/V)

w1 = 2*np.pi*f_axis_1
w2 = 2*np.pi*f_axis_2


#==================RECONSTRUCTED FROM TS PARAMETERS============
loudspeakers = {'full-range1' :'Dayton_CE4895-8',
                'full-range2' :'Dayton_HARB252-8',
                'full-range3' :'SB10PGC21-4',
                'woofer1': 'Peerless_HDSP830860',
                'woofer1_Klippel': 'Peerless_HDSP830860_Klippel',
                'woofer2': 'Dayton_DCS165-4',
                'woofer3': 'Dayton_RS150-4',
                'subwoofer1': 'B&C_15FW76-4',
                'HPtest': 'HPresonant_test'}

loudspeaker = loudspeakers['woofer1_Klippel']

with open(f'Dataset_T&S/{loudspeaker}.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        exec(line)

Z_ts    = Rec + 1j*w1*Lec + Bl**2/(1j*w1*Mms+Rms+1/(1j*w1*Cms))
Hxu_ts  = Bl/((Rec+1j*w2*Lec)*(-(w2**2)*Mms+1j*w2*Rms+1/Cms)+1j*w2*Bl**2)


#=================MAIN PLOT=========================
fig, ax = plt.subplots(1, 2)

ax[0].semilogx(f_axis_1, Z_m, label='measured')
ax[0].semilogx(f_axis_1, Z_f, label='fitted')
ax[0].semilogx(f_axis_1, np.abs(Z_ts), label='reconstructed')
ax[0].set_xlim([10, f_axis_1[-1]])
ax[0].grid(which='both')
ax[0].legend(loc='upper right')
ax[0].set_xlabel('Frequency [Hz]')
ax[0].set_ylabel(r'$|Z(f)|~~[\Omega]$')

ax[1].semilogx(f_axis_2, H_m, label='measured')
ax[1].semilogx(f_axis_2, H_f_noC, label='fitted without creep')
ax[1].semilogx(f_axis_2, H_f, label='fitted')
ax[1].semilogx(f_axis_2, np.abs(Hxu_ts)*1000, label='reconstructed')
ax[1].set_xlim([10, f_axis_2[-1]])
ax[1].grid(which='both')
ax[1].legend(loc='upper right')
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$|H_{XU}(f)|~~[mm/V]$')

plt.show()