#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[set environment]:
# import scipy as sc
from scipy import io
import numpy as np
# import pandas as pd
# import h5py
# import seaborn
# import matplotlib
# import nibabel
import os



# In[set paths etc.]:

workdir='/home/user/Desktop/Projects/bhack_td/'
out_dir=workdir+'results/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)


# In[load analysis parameters]:
fname=out_dir+'bhack_analysis_pars.mat'
tmp=io.loadmat(fname)

tmp['analysis_pars']
analysis_pars=tmp['analysis_pars_cell'][0]

fs                =analysis_pars[0][0][0]  #target sampling frequency for all signals
max_time_s        =analysis_pars[1][0][0]  #maximum time of timecourse to consider for analyses (seconds)
n_chunks          =analysis_pars[2][0][0]  #divide the signal into nchunks "independent" time-courses. 
                                           #GLMs will be cross-validated across these nchunks. For
p_chunk_out       =analysis_pars[3][0][0] #proportion of chunk time points that will be removed
#from the beginning (50%)/end(50%) of each chunk. Discarding some of the timepoints connecting the different chunks increases the
#independence of chunks, and betters the generalization. Also potentially useful if chunks are trials and onset/offset
#effects need to be discarded from analysis

ms_lags           =analysis_pars[4][0] #feature-to-brain-lags (ms; 200 ms = the brain represents
#the feature 200 ms after the waveform reaches the tympanic membrane)

ms_time_folding   =analysis_pars[5][0][0] #this is the folding trick for decreasing the size of the temporal
#distance matrix. Instead of computing the distance between timepoint A and timepoint B, we compute the distance between
#timecourse A and timecourse B, of duration ms_time_fold.

do_zscore         =analysis_pars[6][0][0] #standardize distances (independently for each chunk) before computing any GLM

#analysis parameters in number of samples
ns_max            =analysis_pars[7][0][0] #maximum number of samples of the considered timecourse.
ns_chunk          =analysis_pars[8][0][0] #number of samples in each chunk
ns_chunk_out      =analysis_pars[9][0][0] #make ns_chunk_out even
ns_lags           =analysis_pars[10][0] #feature-to-brain lags in n samples
n_lags            =analysis_pars[11][0][0] #number of considered lags
ns_time_folding   =analysis_pars[12][0][0] #folding window in n samples





# In[load data]:
fname=out_dir+'bhack_task_02_output_temporalfolding.mat'
fname=out_dir+'bhack_task_01_output.mat'
tmp=io.loadmat(fname)


Xs = tmp['Xs'][0]# Xs[0] is acoustics models; Xs[1] is semantic models; Xs[2] is iEEG
Ts = tmp['Ts'][0] #time vectors in ms
Ns=tmp['Ns'][0] #model/sensor names (third dimension of Xs)
fss = tmp['fss'][0]#sampling frequencies




# In[do the temporal folding]:

if ns_time_folding>0:
    print('OK let''s fold the time courses!')
    n_time_foldings=np.floor(Xs[2].shape[0]/ns_time_folding); #number of "foldings"
    
    foldings_idx=np.arange(ns_time_folding*n_time_foldings).astype(int)
    tmp=np.take(Xs[2],foldings_idx,axis=0)
    stmp=tmp.shape
    tmp=np.reshape(tmp, [stmp[0],1,stmp[1],stmp[2],stmp[3]], order='F')
    stmp=tmp.shape
    tmp=np.reshape(tmp,[ns_time_folding.astype(int),n_time_foldings.astype(int),stmp[2],stmp[3],stmp[4]])
    tmp=np.transpose(tmp,[1,0,2,3,4])  #[:,:,:,:,0];
    s=tmp.shape
    tmp=np.reshape(tmp,[s[0],s[1]*s[2],s[3],s[4]],order='F')
    Xs[2]=tmp;
    print('iEEG folded')
    
    
    for acosem in [0,1]:
            for ilag in range(len(ns_lags)):
                ns_thislag=ns_lags[ilag]
                tmp=Xs[acosem]
                tmp=np.take(tmp,foldings_idx+ns_thislag,axis=0)
                stmp=tmp.shape
                tmp=np.reshape(tmp, [stmp[0],1,stmp[1],stmp[2],stmp[3]], order='F')
                stmp=tmp.shape
                tmp=np.reshape(tmp,[ns_time_folding.astype(int),n_time_foldings.astype(int),stmp[2],stmp[3],stmp[4]])
                tmp=np.transpose(tmp,[1,0,2,3,4]);
                s=tmp.shape
                tmp=np.reshape(tmp,[s[0],s[1]*s[2],s[3],s[4]],order='F')
    
                if ilag==0:
                    tmpout=tmp;
                else:
                    tmpout=np.concatenate((tmpout,tmp),axis=2)
            Xs[acosem]=tmpout
    print('Acoustics and semantics folded!')



