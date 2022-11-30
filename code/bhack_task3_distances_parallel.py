# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 10:40:57 2022

@author: marieplgt, chris-zielinski

some multiprocessing tests

"""

import scipy as sc
from scipy import io
import numpy as np
import time
import multiprocessing as mp
import itertools

# initial function with 2 nested for loop
def compute_dist(dat, type='cosine'):
    """Return distance using scipy function
    Parameters
    ----------
    dat : multidimensional array
    type: str
        type of distance: 'cosine' or 'euclidian'
        
    Returns
    -------
    dist: the computed distance
    """
    # get the data shape
    ndim = np.shape(dat)
    n_pair = int(ndim[0]*(ndim[0]-1)/2)
    dist = np.zeros((n_pair, ndim[2], ndim[3]))
    for i in range(ndim[2]):
        for j in range(ndim[3]):
            dist[:, i, j] = sc.spatial.distance.pdist(dat[:,:,i,j], metric=type)
    return dist

def comp_distc(vec, params):
    i = params[0]
    j = params[1]
    d = sc.spatial.distance.pdist(vec[:,:,i,j], metric='euclidean')
    return d   

if __name__ == '__main__':
    # fname = os.path.join('.\results', 'bhack_task_02_output_temporalfolding.mat') 
    fname = 'C:\\home\oprojects\\brainhack22\\bhack_td\\results\\bhack_task_02_output_temporalfolding.mat'
    tmp = io.loadmat(fname)

    Xs = tmp['Xs'][0]
    
    start_time = time.perf_counter()
    dist_eeg = compute_dist(Xs[2], 'euclidean')
    finish_time = time.perf_counter()
    print("Program finished in {} s".format(finish_time-start_time))
    
    start_time = time.perf_counter()
    pool = mp.Pool(mp.cpu_count())
    dat = Xs[2]
    ndim = np.shape(dat)
    comb = list(itertools.product(range(ndim[2]), range(ndim[3])))
    results = [pool.apply(comp_distc, args=(dat, params)) for params in comb]
    finish_time = time.perf_counter()
    print("Program finished in {} s".format(finish_time-start_time))
    pool.close() 
    # Program finished in 0.5110705999977654 s
    # Program finished in 5.315561700001126 s
    # --> it failed...














    
    