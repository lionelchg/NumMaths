import os
import numpy as np

def create_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def make_times(tmin, tend, dt):
    """ Make time vectors from tmind, tend, dt """ 
    times = [tmin]
    while (times[-1] < tend):
        times.append(times[-1] + dt)
    return np.array(times)