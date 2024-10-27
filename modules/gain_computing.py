import numpy as np

def limiter(abs_x, thres):
    return np.minimum(1, thres/abs_x)

def limiter_knee(abs_x, thres, knee):
    return np.minimum(1, 1/(1 + (abs_x/thres)**(1/knee))**knee)

'''
def compressor(abs_x, thres, ratio, knee):
    return np.minimum(1, (thres + (abs_x-thres)/ratio)/abs_x)
'''

def compressor(abs_x, thres, ratio, knee):
    return np.minimum(1, (thres*ration + abs_x - thres)/(ratio*abs_x))


def compressor_knee(abs_x, thres, ratio, knee):
    if 2*(abs_x-thres) < -knee:
        y = abs_x
    elif 2*abs(abs_x-thres) <= knee:
        y = abs_x + (1/ratio - 1)*(x-thres + knee/2)**2/(2*knee)
    else:
        y = thres + (abs_x-thres)/ratio
    return y/abs_x