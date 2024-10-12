import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from lib import multibandEllipticFilters, multibandChebyshev1Filters, multibandChebyshev2Filters, multibandBesselFilters


Fs = 48000
f  = np.geomspace(20, Fs/2, 1000)
w  = 2*np.pi*f

#create a 5-band equalizer with remez algorithm

# Multi-band processing
fc = [40, 90, 200]     # cutoff frequencies

attack_time  = 0.00005      # Attack time in seconds
release_time = 0.07         # Release time in seconds
averaging_time = 0.05       # Averaging time in seconds

sos_filters = multibandBesselFilters(4, fc, Fs, norm='phase')
# sos_filters = multibandChebyshev1Filters(4, fc, Fs, 1)
# sos_filters = multibandChebyshev2Filters(4, fc, Fs, 20)
#sos_filters = multibandEllipticFilters(4, fc, Fs, 1, 60)

H_filter = np.zeros((len(sos_filters), len(f)), dtype=complex)
# Group delay
G_filter = np.zeros((len(sos_filters), len(f)))

for i in range(len(sos_filters)):
    _, H_filter[i] = sig.sosfreqz(sos_filters[i], worN=f, fs=Fs)
    
    # Compute group delay by differentiating the unwrapped phase with respect to frequency
    phase = np.unwrap(np.angle(H_filter[i]))
    G_filter[i] = -np.gradient(phase, w)  # Negative derivative of phase w.r.t frequency

# Compute total response by combining filters
H_tot = H_filter[0] - H_filter[1] + H_filter[2] - H_filter[3]

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].set_title('Modulus and Group Delay of the Filters')

for i in range(len(sos_filters)):
    ax[0].semilogx(f, 20 * np.log10(np.abs(H_filter[i])), label=f'Band {i+1}')
    ax[1].semilogx(f, G_filter[i], label=f'Band {i+1}')

# Plot total filter response
#ax[0].semilogx(f, 20 * np.log10(np.abs(H_tot)), label='Total')
#ax[1].semilogx(f, np.unwrap(np.angle(H_tot)), label='Total')

ax[0].set(xlabel='Frequency', ylabel='Magnitude [dB]')
ax[1].set(xlabel='Frequency', ylabel='Group Delay [s]')

ax[0].set(xlim=(f[0], Fs/2), ylim=(-40, 10))

ax[0].legend()
ax[1].legend()

plt.show()
