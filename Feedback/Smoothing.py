import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 19})
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

duration = 0.3 #[s]
Fs = 48000
t = np.arange(0, duration, 1/Fs)

# create a step function

t_start = 0.025
t_end = 0.1

Cms_comp = np.ones(len(t))
Cms_target = np.ones(len(t))
Cms_target[int(t_start*Fs):int(t_end*Fs)] = 0

attack = 0.03
release = 0.1

attack_coeff = 1 - np.exp(-2.2 / (attack * Fs))
release_coeff = 1 - np.exp(-2.2 / (release * Fs))

tau_attack = attack/2.2
tau_release = release/2.2

for n in range(1, len(Cms_target)):

    if Cms_comp[n] > Cms_target[n]:
        k = attack_coeff
    else:
        k = release_coeff

    Cms_comp[n] = (1 - k) * Cms_comp[n-1] + k * Cms_target[n]


fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')
ax.plot(t*1e3, Cms_target, color='0.2', label=r'$C_{ms}^{(target)}[n]$')
ax.plot(t*1e3, Cms_comp, color='0.0',label=r'$C_{ms}^{(comp)}[n]$', linestyle='--',dashes=(7, 5))

#add vertical lines correspondind to the 10% and 90% values of the step response
t_fall_start = t_start - np.log(0.9)*tau_attack
t_fall_end = t_start - np.log(0.1)*tau_attack

t_rise_start = t_end - np.log(0.9)*tau_release
t_rise_end = t_end - np.log(0.1)*tau_release

ax.axvline(t_fall_start*1e3, color='0.5', linestyle=':')
ax.axvline(t_fall_end*1e3, color='0.5', linestyle=':')
ax.axvline(t_rise_start*1e3, color='0.5', linestyle=':')
ax.axvline(t_rise_end*1e3, color='0.5', linestyle=':')

#add horizontal arrow to between t_fall_start and t_fall_end and for t_rise_start and t_rise_end with 't_a' and 't_r' labels
#above the middle of the arrows

ax.annotate('', xy=(t_fall_end*1e3, 0.5), xytext=(t_fall_start*1e3, 0.5), arrowprops=dict(arrowstyle='<->'))
ax.annotate(r'$t_a$', xy=(t_fall_start*1e3 + (t_fall_end*1e3 - t_fall_start*1e3)/2 - 2, 0.53))

ax.annotate('', xy=(t_rise_end*1e3, 0.5), xytext=(t_rise_start*1e3, 0.5), arrowprops=dict(arrowstyle='<->'))
ax.annotate(r'$t_r$', xy=(t_rise_start*1e3 + (t_rise_end*1e3 - t_rise_start*1e3)/2 - 2, 0.53))

ax.set(xlabel='Time [ms]', ylabel=r'$C_{ms}$ [m/N]')
ax.legend(loc='lower right')
ax.grid()
plt.show()