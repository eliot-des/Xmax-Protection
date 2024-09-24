import numpy as np
import scipy.signal as sig

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