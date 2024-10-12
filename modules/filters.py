import numpy as np
import scipy.signal as sig


class EqualizerFilter:

    def __init__(self, f_type, fc, Q, dBgain=0, fs=48e3):
        # f_type ("LP","HP")
        self.b, self.a = self.get_filter_coefficients(f_type, fc, Q, dBgain, fs)

    def get_filter_coefficients(self, f_type, Fc, Q, dBgain, fs):

        wc = 2*np.pi*Fc/fs
        alpha = np.sin(wc)/(2*Q)
        A = 10**(dBgain/40)


        if f_type == "LP":

            b0 = (1 - np.cos(wc))/2
            b1 = 1 - np.cos(wc)
            b2 = (1 - np.cos(wc))/2
            a0 = 1 + alpha
            a1 = -2*np.cos(wc)
            a2 = 1 - alpha

        elif f_type == "HP":

            b0 = (1 + np.cos(wc))/2
            b1 = -(1 + np.cos(wc))
            b2 = (1 + np.cos(wc))/2
            a0 = 1 + alpha
            a1 = -2*np.cos(wc)
            a2 = 1 - alpha

        elif f_type == "BP":

            b0 = alpha
            b1 = 0
            b2 = -alpha
            a0 = 1 + alpha
            a1 = -2*np.cos(wc)
            a2 = 1 - alpha

        elif f_type == "peak":

            b0 = 1 + alpha*A
            b1 = -2*np.cos(wc)
            b2 = 1 - alpha*A
            a0 = 1 + alpha/A
            a1 = -2*np.cos(wc)
            a2 = 1 - alpha/A

        elif f_type == "LS":

            b0 = A*((A+1) - (A-1)*np.cos(wc) + 2*np.sqrt(A)*alpha)
            b1 = 2*A*((A-1) - (A+1)*np.cos(wc))
            b2 = A*((A+1) - (A-1)*np.cos(wc) - 2*np.sqrt(A)*alpha)
            a0 = (A+1) + (A-1)*np.cos(wc) + 2*np.sqrt(A)*alpha
            a1 = -2*((A-1) + (A+1)*np.cos(wc))
            a2 = (A+1) + (A-1)*np.cos(wc) - 2*np.sqrt(A)*alpha

        elif f_type == "HS":

            b0 = A*((A+1) + (A-1)*np.cos(wc) + 2*np.sqrt(A)*alpha)
            b1 = -2*A*((A-1) + (A+1)*np.cos(wc))
            b2 = A*((A+1) + (A-1)*np.cos(wc) - 2*np.sqrt(A)*alpha)
            a0 = (A+1) - (A-1)*np.cos(wc) + 2*np.sqrt(A)*alpha
            a1 = 2*((A-1) - (A+1)*np.cos(wc))
            a2 = (A+1) - (A-1)*np.cos(wc) - 2*np.sqrt(A)*alpha

        elif f_type == "notch":

            b0 = 1
            b1 = -2*np.cos(wc)
            b2 = 1
            a0 = 1 + alpha
            a1 = -2*np.cos(wc)
            a2 = 1 - alpha

        elif f_type == "AP":

            b0 = 1 - alpha
            b1 = -2*np.cos(wc)
            b2 = 1 + alpha
            a0 = 1 + alpha
            a1 = -2*np.cos(wc)
            a2 = 1 - alpha

        else:
            raise ValueError('The filtre type {} is incorrect'.format(f_type))

        # normalize filter coefficients so that a0 = 1
        b = [b0/a0, b1/a0, b2/a0]
        a = [1, a1/a0, a2/a0]

        return b, a


class PeakFilterTD2:
    def __init__(self, Fc, Q, dBgain=0, fs=48e3):

        self.fs = fs
        self.get_filter_coefficients(Fc, Q, dBgain)
        
        # memory
        
        self.d0 = 0
        self.d1 = 0
        self.d2 = 0
        self.d3 = 0

    def get_filter_coefficients(self, Fc, Q, dBgain):

        wc = 2*np.pi*Fc/self.fs
        alpha = np.sin(wc)/(2*Q)
        A = 10**(dBgain/40)

        b0 = 1 + alpha*A
        b1 = -2*np.cos(wc)
        b2 = 1 - alpha*A
        a0 = 1 + alpha/A
        a1 = -2*np.cos(wc)
        a2 = 1 - alpha/A

        # normalize filter coefficients so that a0 = 1
        self.b = [b0/a0, b1/a0, b2/a0]
        self.a = [1,     a1/a0, a2/a0]

    def filter(self, x, Fc, Q, dBgain):
        #recompute filter coefficients
        self.get_filter_coefficients(Fc, Q, dBgain)

        y = self.b[0]*x + self.b[1]*self.d0 + self.b[2]*self.d1 - self.a[1]*self.d2 - self.a[2]*self.d3

        #update memory
        self.d1 = self.d0 
        self.d0 = x
        self.d3 = self.d2
        self.d2 = y

        return y

class HighPassFilterTD2:
    def __init__(self, Fc, Q, fs=48e3):
        self.fs = fs
        self.z1 = 0
        self.z2 = 0
        self.get_filter_coefficients(Fc, Q)

    def get_filter_coefficients(self, Fc, Q):
        wc = 2*np.pi*Fc/self.fs
        alpha = np.sin(wc)/(2*Q)
        b0 = (1 + np.cos(wc))/2
        b1 = -(1 + np.cos(wc))
        b2 = (1 + np.cos(wc))/2
        a0 = 1 + alpha
        a1 = -2 * np.cos(wc)
        a2 = 1 - alpha

        # Normalize filter coefficients so that a0 = 1
        self.b = [b0/a0, b1/a0, b2/a0]
        self.a = [1, a1/a0, a2/a0]

    def filter(self, x, Fc, Q):
        # Recompute filter coefficients
        self.get_filter_coefficients(Fc, Q)

        # Direct Form II Transposed implementation
        y = self.b[0]*x + self.z1
        self.z1 = self.b[1]*x - self.a[1]*y + self.z2
        self.z2 = self.b[2]*x - self.a[2]*y

        return y