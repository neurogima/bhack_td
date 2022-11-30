## Do the imports
import os
import numpy as np
import scipy.io as io

## set paths and analysis parameters: run for all tasks
rootmain='/home/jnfn/Work/Projects/Code/Sound_Brain_Sem/bhack_td/'
out_dir='results/'
fname=rootmain+out_dir+'bhack_analysis_pars.mat'
tmp=io.loadmat(fname)
tmp=['analysis_pars']
analysis_pars=tmp['analysis_pars_cell'][0]

analysis_pars=tmp['analysis_pars_cell'][0]

fs                =analysis_pars[0][0][0]  #target sampling frequency for all signals
max_time_s        =analysis_pars[1][0][0]  #maximum time of timecourse to consider for analyses (seconds)
n_chunks          =analysis_pars[2][0][0]  #divide the signal into nchunks "independent" time-courses. 
                                           #GLMs will be cross-validated across these nchunks. For
p_chunk_out       =analysis_pars[3][0][0] #proportion of chunk time points that will be removed
#from the beginning (50%)/end(50%) of each chunk. Discarding some of the timepoints connecting the different chunks increases the
#independence of chunks, and betters the generalization. Also potentially useful if chunks are trials and onset/offset
#effects need to be discarded from analysis

ms_lags           =analysis_pars[4][0][0] #feature-to-brain-lags (ms; 200 ms = the brain represents
#the feature 200 ms after the waveform reaches the tympanic membrane)

ms_time_folding   =analysis_pars[5][0][0] #this is the folding trick for decreasing the size of the temporal
#distance matrix. Instead of computing the distance between timepoint A and timepoint B, we compute the distance between
#timecourse A and timecourse B, of duration ms_time_fold.

do_zscore         =analysis_pars[6][0][0] #standardize distances (independently for each chunk) before computing any GLM

#analysis parameters in number of samples
ns_max            =analysis_pars[7][0][0] #maximum number of samples of the considered timecourse.
ns_chunk          =analysis_pars[8][0][0] #number of samples in each chunk
ns_chunk_out      =analysis_pars[9][0][0] #make ns_chunk_out even
ns_lags           =analysis_pars[10][0][0] #feature-to-brain lags in n samples
n_lags            =analysis_pars[11][0][0] #number of considered lags
ns_time_folding   =analysis_pars[12][0][0] #folding window in n samples

fnout_tasks = ['bhack_task_01_output.mat', 'bhack_task_02_output_temporalfolding.mat', 'bhack_task_03_output_distances.mat'
    ...: , 'bhack_task_04_output_GLMs.mat', 'bhack_task_05_output_plots.pdf']
    
fnout_tasks = [os.path.join(rootmain, out_dir, name) for name in fnout_tasks]

##In[permute EEG data]:

nobj=Y.shape[2]
nperms=25;
s = np.random.default_rng().uniform(-1,0,[nobj,nperms])
perms=np.argsort(s,axis=0)

for i in range(nperms):
    permY=Y[:,:,perms[:,i],: ]


# -*- coding: utf-8 -*-

