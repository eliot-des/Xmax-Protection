import numpy as np
import matplotlib.pyplot as plt

# Develop a state-space model for a loudspeaker

# Thiele-small parameters
# Data from Peerless HDS-P830860 Datasheet

fs = 72
Rec = 6.4
Lec = 0.278e-3

Qms = 2.08
Qes = 0.725
Qts = 0.54

Mms = 8.88e-3
Cms = 560e-6
Rms = 1/(2*np.pi*fs*Cms*Qms)

Bl  = 5.74

Vas = 6.36e-3
Sd  = 89.9e-4


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
t = np.arange(0, 1, Ts)
u = np.sin(2*np.pi*f*t)


x = np.zeros((3, len(t)))

for i in range(0, len(t)-1):
    x[:,i+1] = Ad @ x[:,i] + (Bd * u[i])[:,0]


fig, ax = plt.subplots()

ax.plot(t, x[0], label='x1 - current')
ax.plot(t, 1e3*x[1], label='x2 - displacement')
ax.plot(t, x[2], label='x3 - velocity')

ax.legend()
plt.show()



# Hxu - displacement transfer function

f = np.arange(20, 20000)
w = 2*np.pi*f

Hxu = Bl/((Rec+1j*w*Lec)*(-(w**2)*Mms+1j*w*Rms+1/Cms)+1j*w*Bl**2)


fig, ax = plt.subplots()

ax.semilogx(f, np.abs(Hxu)*1e3)
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('|X/U| [mm/V]')
ax.grid(which='both')

plt.show()