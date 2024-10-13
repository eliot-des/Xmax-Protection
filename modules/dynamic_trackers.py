import numpy as np
import scipy.signal as sig

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

def rms_follower_sbs(x, x_prev, averaging_time, Fs):

    averaging_coeff = 1 - np.exp(-2.2 / (averaging_time * Fs))

    #squared RMS value
    x_rms_square = (1 - averaging_coeff) * x_prev + averaging_coeff * x**2
    # x_rms_output = np.sqrt(x2_rms)
    
    return x_rms_square


def gain_factor_smoothing(x, attack_time, release_time, Fs):
    y = np.zeros_like(x)

    attack_coeff  = 1 - np.exp(-2.2 / (attack_time * Fs))
    release_coeff = 1 - np.exp(-2.2 / (release_time * Fs))

    for n in range(1, len(x)):
        if x[n] > y[n-1]:
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