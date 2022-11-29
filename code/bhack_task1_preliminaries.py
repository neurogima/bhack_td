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

dat_dir=workdir+'dat/'

fn_ieeg=dat_dir+'brain.mat' #ieeg data
fn_aco=dat_dir+'aco.mat'   #all acoustic models
fn_sem=dat_dir+'sem.mat'   #semantic models
dat_fns=[fn_ieeg,fn_aco,fn_sem]
dat_nams=['eeg','aco','sem']
var_nams=['X','fs','names','time']

# The content of each dataset is four variables
# X = the signal of size [ntimepoints nembeddingdimensions n_sensors(iEEG)_OR_nmodels(aco/sem)]
# fs = the sampling frequency in Hz
# names = the names for the sensors/models (third dimension of X)
# time = time vector in ms
    
for i in range(len(dat_fns)):
    tmp=io.loadmat(dat_fns[i])
    print('---------------')
    print(dat_fns[i])
    for j in range(len(var_nams)): 
        # #%%%% loop through the variables contained in this .mat dataset
        varnam=var_nams[j]
        varnew=varnam+'_'+dat_nams[i] #%new name for the variable allocated to workspace - appends datasat name (eeg/aco/sem)
        exec(varnew+"=tmp['"+varnam+"']")


for i in range(len(names_aco)):
    names_aco[i]=names_aco[i][0][0]
for i in range(len(names_sem)):
    names_sem[i]=names_sem[i][0][0]
for i in range(len(names_eeg)):
    names_eeg[i]=names_eeg[i][0][0]




Xs = [X_aco,X_sem,X_eeg] # %datasets
Ts = [np.round(time_aco,decimals=0),np.round(time_sem,decimals=0),np.round(time_eeg,decimals=0)] #time vectors in ms
#Ts = celfun(@round,Ts); %round to the nearest millisecond, eliminates very small timing vector gremlins
Ns = [names_aco,names_sem,names_eeg] # %model/sensor names (third dimension of Xs)
fss = [fs_aco,fs_sem,fs_eeg] # %sampling frequencies


## In[deal with infs in acoustics and nans in semantics]:

Xs[0][np.isinf(Xs[0])]=0 #replace inf values in acoustics with zeroes;


# replace time samples without semantic embedding with average embedding 
# of other time samples. This is a temporary hack, should be done independently
# for each cross-validation fold (i.e., temporal chunk)
for i in [0,1]:
    tmp=Xs[1][:,:,i]
    tmpnanmean=np.nanmean(tmp,axis=0,keepdims=True)
    idxnan=np.isnan(tmp[:,1])
    tmp[idxnan,:]=tmpnanmean
    Xs[1][:,:,i]=tmp

    


## In[downsample time courses]:

# let's do a very crude downsampling of data to the new sampling frequency
for i in range(len(Xs)):
    tmpfs = fss[i] # %the sampling frequency of this signal in Xs{i}
    tmpdown = tmpfs/fs #; %down-sampling factor: take one sample every tmpdown samples
    tmpdown = tmpdown.astype(int)[0][0]
    
    #if rem(tmpdown,1)~=0 || tmpdown<1 %if not integer or if target fs > tmpfs (upsampling)
    #     error('the ratio of the old fs to the new fs must be an integer; target fs must be < old fs')
    #end
    
    #apply the downsampling;
    Xs[i]=Xs[i][::tmpdown,:,:]
    Ts[i]=Ts[i][::tmpdown]

print([x.shape for x in Xs])





# In[add differential to acoustics features]:
tmpacodiff=np.diff(Xs[0],n=1,prepend=0)
Xs[0]=np.concatenate((Xs[0],tmpacodiff),axis=2)

newnames=names_aco.copy()
for i in range(len(newnames)):
    newnames[i]=newnames[i]+'Diff'
Ns[0]=np.concatenate((Ns[0],newnames),axis=0)


print([x.shape for x in Xs])


# In[cut all timecourses to same duration]:


for i in range(len(Xs)):
    
    thisns=np.min([ns_max,Ts[i].shape[0]])
    Xs[i]=Xs[i][0:thisns,:,:]
    Ts[i]=Ts[i][0:thisns,:]

if np.sum(np.equal(Ts[0],Ts[1]))==thisns:
    if np.sum(np.equal(Ts[0],Ts[2]))==thisns:
        print('All good with time vectors! let''s continue')

## In[create the temporal chunks]:


chunks_idx=np.arange(Ts[0].shape[0]).astype(int)
chunks_idx=np.reshape(chunks_idx, [ns_chunk,n_chunks], order='F')

for i in range(len(Xs)):
    for j in range(n_chunks):
        tmp=Xs[i][chunks_idx[:,j],:,:]
        tmp=np.expand_dims(tmp,3) #add one last singletone dimension
        if j==0:
            tmpout=tmp
        else:
            tmpout=np.concatenate((tmpout,tmp),axis=-1)
    Xs[i]=tmpout

print([x.shape for x in Xs])



# In[chop out more samples for independence of data from different temporal chunks etc.]:

#let's chop out ns_chunk_out/2 from beginning and end. This aids the
#independence of data in different trials.
inidx=np.arange(ns_chunk_out/2+1,Xs[0].shape[0]-ns_chunk_out/2+1,1,dtype=int)
for i in range(len(Xs)):
    Xs[i]=Xs[i][inidx,:,:,:]


print([x.shape for x in Xs])

#let's remove max(ns_lags) from end of eeg (Xs(3)) to account for the
#maximum feature-to-brain lag
inidx=np.arange(np.max(ns_lags),Xs[2].shape[0],1,dtype=int)

Xs[2]=Xs[2][inidx,:,:,:];

print([x.shape for x in Xs])






