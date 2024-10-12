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





def normalize(data):
    #normalize data depending on the bit depth
    if data.dtype == np.int16:
        data = data/2**15
    elif data.dtype == np.int32:
        data = data/2**31
    elif data.dtype == np.float32:
        data = data
    else: 
        raise ValueError('Bit depth not supported')

    data = data/np.max(np.abs(data))
    
    return data


#function to create a band-pass filter using two butterworth filters

def linkwitzBandPass(order, fc1, fc2, fs):
    if order % 2 != 0:
        raise ValueError('Order must be even')
    else:
        sos_lp = sig.butter(order, fc2, 'low', fs=fs, output='sos')
        sos_hp = sig.butter(order, fc1, 'high', fs=fs, output='sos')

        sos_bp = np.concatenate((sos_lp, sos_hp))

        # Squaring the transfer function
        sos_bp = np.concatenate((sos_bp, sos_bp))

    return sos_bp
#function to create a low-pass Linkwitz-Riley filter
def linkwitzLowPass(order, fc, fs):
    if order%2 != 0:
        raise ValueError('Order must be even')
    else:
        sos_lp = sig.butter(order, fc, 'low', fs = fs, output='sos')
        sos_lp = np.concatenate((sos_lp, sos_lp))
    return sos_lp


#function to create a high-pass Linkwitz-Riley filter
def linkwitzHighPass(order, fc, fs):
    if order%2 != 0:
        raise ValueError('Order must be even')
    else:
        sos_hp = sig.butter(order, fc, 'high', fs = fs, output='sos')
        sos_hp = np.concatenate((sos_hp, sos_hp)) 
    return sos_hp

#function to return a list of sos sections for multi-band processing with Linkwitz-Riley filters
def multibandLinkwitzFilters(order, fc, fs):
    
    Nbrfilters = len(fc) +1

    sos = [None] * Nbrfilters

    #create the low and high pass filters
    sos[0] = linkwitzLowPass(order, fc[0], fs)
    sos[-1] = linkwitzHighPass(order, fc[-1], fs)

    #create the band-pass filters
    for i in range(1, Nbrfilters-1):
        sos[i] = linkwitzBandPass(order, fc[i-1], fc[i], fs)
    
    return sos


##function to return a list of sos sections for multi-band processing with elliptic filters
def multibandEllipticFilters(order, fc, fs, rp, rs):
        
        Nbrfilters = len(fc) +1
    
        sos = [None] * Nbrfilters
    
        #create the low and high pass filters
        sos[0] = sig.ellip(order, rp, rs, fc[0], 'low', fs = fs, output='sos')
        sos[-1] = sig.ellip(order, rp, rs, fc[-1], 'high', fs = fs, output='sos')
    
        #create the band-pass filters
        for i in range(1, Nbrfilters-1):
            sos[i] = sig.ellip(order, rp, rs, (fc[i-1], fc[i]), 'bandpass', fs = fs, output='sos')

        return sos


##function to return a list of sos sections for multi-band processing with butterworth filters
def multibandButterworthFilters(order, fc, fs):
    
    Nbrfilters = len(fc) +1

    sos = [None] * Nbrfilters

    #create the low and high pass filters
    sos[0] = sig.butter(order, fc[0], 'low', fs = fs, output='sos')
    sos[-1] = sig.butter(order, fc[-1], 'high', fs = fs, output='sos')

    #create the band-pass filters
    for i in range(1, Nbrfilters-1):
        sos[i] = sig.butter(order, (fc[i-1], fc[i]), 'bandpass', fs = fs, output='sos')

    return sos

#function to return a list of sos sections for multi-band processing with Chebyshev filters
def multibandChebyshev1Filters(order, fc, fs, rp):
    
    Nbrfilters = len(fc) +1

    sos = [None] * Nbrfilters

    #create the low and high pass filters
    sos[0] = sig.cheby1(order, rp, fc[0], 'low', fs = fs, output='sos')
    sos[-1] = sig.cheby1(order, rp, fc[-1], 'high', fs = fs, output='sos')

    #create the band-pass filters
    for i in range(1, Nbrfilters-1):
        sos[i] = sig.cheby1(order, rp, (fc[i-1], fc[i]), 'bandpass', fs = fs, output='sos')
        
    return sos

#function to return a list of sos sections for multi-band processing with Chebyshev type 2 filters
def multibandChebyshev2Filters(order, fc, fs, rs):

    Nbrfilters = len(fc) +1

    sos = [None] * Nbrfilters

    #create the low and high pass filters
    sos[0] = sig.cheby2(order, rs, fc[0], 'low', fs = fs, output='sos')
    sos[-1] = sig.cheby2(order, rs, fc[-1], 'high', fs = fs, output='sos')

    #create the band-pass filters
    for i in range(1, Nbrfilters-1):
        sos[i] = sig.cheby2(order, rs, (fc[i-1], fc[i]), 'bandpass', fs = fs, output='sos')
        
    return sos

#function to return a list of sos sections for multi-band processing with bessel filters
def multibandBesselFilters(order, fc, fs, norm='phase'):
        
        Nbrfilters = len(fc) +1
    
        sos = [None] * Nbrfilters
    
        #create the low and high pass filters
        sos[0] = sig.bessel(order+1, fc[0], 'low', norm=norm, fs = fs, output='sos')
        sos[-1] = sig.bessel(order+1, fc[-1], 'high', norm=norm, fs = fs, output='sos')
    
        #create the band-pass filters
        for i in range(1, Nbrfilters-1):
            sos[i] = sig.bessel(order, (fc[i-1], fc[i]), 'bandpass', fs = fs, output='sos')
            
        return sos



def dynamic_peak_follower1(x, attack_time, release_time, sample_rate):
    x_peak = 0.0
    x_peak_output = np.zeros_like(x)

    # Calculate attack and release coefficients based on provided times
    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * sample_rate))
    release_coeff = 1 - np.exp(-2.2 / (release_time * sample_rate))

    for n in range(len(x)):
        abs_x = np.abs(x[n])
        
        delta = abs_x - x_peak
        
        if delta > 0:
            x_peak += attack_coeff * delta
        else:
            x_peak += release_coeff * delta
        
        # Store the result in the output array
        x_peak_output[n] = x_peak

    return x_peak_output

def dynamic_peak_follower2(x, attack_time, release_time, Fs):
    
    x_peak = 0  
    x_peak_output = np.zeros_like(x)
    
    # Attack and release coefficients based on time constants
    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
    release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

    # Process each sample
    for n in range(len(x)):

        abs_x = np.abs(x[n])
        
        if abs_x > x_peak:
            x_peak = (1 - attack_coeff) * x_peak + attack_coeff * abs_x
        else:
            x_peak = (1 - release_coeff) * x_peak
        
                
        # Store the result
        x_peak_output[n] = x_peak
    
    return x_peak_output


#sample by sample peak follower with attack and release time constants
def dynamic_peak_follower_sbs(x, x_peak, attack_time, release_time, Fs):
        
        # Attack and release coefficients based on time constants
        attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
        release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))
    
        # Process the sample
        abs_x = np.abs(x)
        
        if abs_x > x_peak:
            x_peak = (1 - attack_coeff) * x_peak + attack_coeff * abs_x
        else:
            x_peak = (1 - release_coeff) * x_peak
    
        return x_peak

def envelope_follower_klippel(x, fr, Fs):

    y = np.zeros_like(x)

    wr = 2*np.pi*fr

    for n in range(1, len(x)):
        dx = (x[n] - x[n-1]) * Fs

        y[n] = np.sqrt(x[n]**2 + (dx/wr)**2)

    return y

#RMS Measurement
def dynamic_rms_follower(x, averaging_time, Fs):
    
    x2_rms = 0  
    x_rms_output = np.zeros_like(x)

    averaging_coeff = 1 - np.exp(-2.2 / (averaging_time * Fs))

    for n in range(len(x)):
        #squared RMS value
        x2_rms = (1 - averaging_coeff) * x2_rms + averaging_coeff * x[n]**2
        x_rms_output[n] = np.sqrt(x2_rms)
    
    return x_rms_output


def gain_factor_smoothing(x, attack_time, release_time, Fs):
    y = np.zeros_like(x)

    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
    release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

    for n in range(1, len(x)):
        if x[n] > x[n-1]:
            k = attack_coeff
        else:
            k = release_coeff

        y[n] = (1 - k) * y[n-1] + k * x[n]

    return y


def gain_factor_smoothing_sbs(x, x_initial, attack_time, release_time, Fs):

    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
    release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

    if x > x_initial:
        k = attack_coeff
    else:
        k = release_coeff
    
    x = (1 - k) * x_initial + k * x

    return x


def gain_factor_smoothing_sbs_bis(x, x_initial, attack_time, release_time, Fs):

    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
    release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

    if x < x_initial:    # change of the sign compared to gain_factor_smoothing_sbs
        k = attack_coeff
    else:
        k = release_coeff
    
    x = (1 - k) * x_initial + k * x

    return x