import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqs, freqz

# Develop a state-space model for a loudspeaker

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

# State-space model
# Create the Matrices A, B and vectors x, xdot

x    = np.array([[0], [0], [0]])
xdot = np.array([[0], [0], [0]])

A =   np.array([[-Rec/Lec,          0,    -Bl/Lec], 
                [0,                 0,          1],
                [Bl/Mms, -1/(Cms*Mms),   -Rms/Mms]])

B =  np.array([[1/Lec], [0], [0]])

I = np.eye(3)

# Solve system for x[n+1] = Ad*x[n] + Bd*u[n]
# Bilinear transform

Fs = 44100
Ts = 1/Fs


Mtemp = np.linalg.inv(I - (Ts/2)*A)

Ad = Mtemp @ (I + (Ts/2)*A)
Bd = (Mtemp @ B) * Ts

print(Bd.shape)


f = 200
t = np.arange(0, 0.02, Ts)
u = np.sin(2*np.pi*f*t)


x = np.zeros((3, len(t)))

for i in range(0, len(t)-1):
    x[:,i+1] = Ad @ x[:,i] + (Bd * u[i])[:,0]


fig, ax = plt.subplots()

ax.plot(t, x[0], label='x1 - current')
ax.plot(t, 1e3*x[1], label='x2 - displacement')
ax.plot(t, x[2], label='x3 - velocity')

ax.legend()





# Hxu - displacement transfer function

f = np.geomspace(20, 20000, 500)
w = 2*np.pi*f
s = 1j*w

C = np.array([[0, 1, 0]])
D = np.array([[0]])


#analog transfer function
H = np.array([C @ np.linalg.inv((s_*I - A)) @ B + D for s_ in s])[:,0,0]

Hxu = Bl/((Rec+1j*w*Lec)*(-(w**2)*Mms+1j*w*Rms+1/Cms)+1j*w*Bl**2)
Ze  = Rec+1j*w*Lec + Bl**2/(Rms + 1j*w*Mms + 1/(1j*w*Cms))
# Filter coefficient from S to z domain

b = np.array([0, 0, 0, Bl])
a = np.array([Lec*Mms, Rec*Mms+Lec*Rms, Rec*Rms+Lec/Cms+Bl**2, Rec/Cms])

#normalize filter
a0 = a[0]
a = a/a0
b = b/a0


_, h = freqs(b, a, worN=w)

# Filter coefficient from S to z domain analytically

B_LF = np.array([0, 0, Bl/Rec])
A_LF = np.array([Mms, Rms+Bl**2/Rec, 1/Cms])

b_LF = np.zeros(3)
b_LF[0] = B_LF[2]/Fs**2
b_LF[1] = 2*b_LF[0]
b_LF[2] = b_LF[0]

a_LF = np.zeros(3)
a_LF[0] = A_LF[2]/Fs**2 + 2*A_LF[1]/Fs + 4*A_LF[0]
a_LF[1] = 2*A_LF[2]/Fs**2 - 8*A_LF[0]
a_LF[2] = A_LF[2]/Fs**2 + 4*A_LF[0] - 2*A_LF[1]/Fs

_, H_LF = freqs(B_LF, A_LF, worN=w)
_, h_LF = freqz(b_LF, a_LF, worN=f, fs=Fs)


fig, ax = plt.subplots()
ax.semilogx(f, np.abs(H)*1e3, label='from state-space')
ax.semilogx(f, np.abs(Hxu)*1e3, 'r--', label='theoritical')
ax.semilogx(f, np.abs(h)*1e3, '-.m', label='analog filter')
ax.semilogx(f, np.abs(H_LF)*1e3, '-.k', label='analog filter LF approx')
ax.semilogx(f, np.abs(h_LF)*1e3, ':b', label='analytically digital\n filter LF approx')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('|X/U| [mm/V]')
ax.grid(which='both')
ax.legend(loc='upper right')



fig, ax = plt.subplots()

axt = ax.twinx()

ax.semilogx(f, np.abs(Ze))
axt.semilogx(f, np.angle(Ze), 'r')
ax.set(xlabel='Frequency [Hz]', ylabel='Impedance [Ohm]')
axt.set(ylabel='Phase [rad]')
#axt.tick_params(axis='y')
ax.grid(which='both')
plt.show()
