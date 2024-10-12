import numpy as np

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