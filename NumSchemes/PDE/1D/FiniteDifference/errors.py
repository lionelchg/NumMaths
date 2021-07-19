import numpy as np

def L1error(y_th, y_num, ncx):
    """ Compute the L1 error """
    return np.sum(abs(y_th - y_num)) / ncx

def L2error(y_th, y_num, ncx):
    """ Compute the L2 error """
    return np.sqrt(np.sum((y_th - y_num)**2)) / ncx

def Linferror(y_th, y_num, ncx):
    """ Compute the Linf error """
    return np.max(abs(y_th - y_num))

