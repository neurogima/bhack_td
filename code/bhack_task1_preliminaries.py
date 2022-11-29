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




Xs = [X_aco,X_sem,X_eeg] # %datasets: this could be a python list?
Ts = [np.round(time_aco,decimals=0),np.round(time_sem,decimals=0),np.round(time_eeg,decimals=0)] #time vectors in ms
#Ts = celfun(@round,Ts); %round to the nearest millisecond, eliminates very small timing vector gremlins
Ns = [names_aco,names_sem,names_eeg] # %model/sensor names (third dimension of Xs)
fss = [fs_aco,fs_sem,fs_eeg] # %sampling frequencies



#Xs{1}(isinf(Xs{1}))=0; %replace inf values in acoustics with zeroes;
    
# replace time samples without embedding with average embedding of other
# time samples. This is a temporary hack.
#tmpisnan=sum(isnan(X_sem),2)>0;
#    for i=1:2 %give average semantic embedding for time samples without average embedding.
#        %this should be done independently for each chunk
#        tmpnanmean=nanmean(X_sem(:,:,i),1);
#        tmpisnan=sum(isnan(X_sem(:,:,i)),2)>0;
#        Xs{2}(tmpisnan,:,i)=ones(sum(tmpisnan),1)*tmpnanmean;
#    end



# In[]:

# let's do a very crude downsampling of data to the new sampling frequency
for i in range(len(Xs)):
    tmpfs = fss[i] # %the sampling frequency of this signal in Xs{i}
    tmpdown = tmpfs/fs #; %down-sampling factor: take one sample every tmpdown samples
    tmpdown = tmpdown.astype(int)[0][0]
    #if rem(tmpdown,1)~=0 || tmpdown<1 %if not integer or if target fs > tmpfs (upsampling)
    #     error('the ratio of the old fs to the new fs must be an integer; target fs must be < old fs')
    #end
    
    #apply the downsampling;
    #Xs{i} = Xs{i}(1:tmpdown:end,:,:,:); %referencing also 4th dimension for Paul's data
    #Ts{i} = Ts{i}(1:tmpdown:end);
    print(tmpdown)
    Xs[i]=Xs[i][::tmpdown,:,:]
    Ts[i]=Ts[i][::tmpdown]


# In[add differential to acoustics features]:
tmpacodiff=Xs[0];

tmpzeroes=np.zeros(tmpaco[0:1,:,:].shape)
tmpacodiff=np.concatenate((tmpzeroes,tmpaco),axis=0)
tmpacodiff=tmpacodiff[1:-1,:,:]-tmpacodiff[0:-2,:,:]
Xs[0]=np.concatenate((Xs[0],tmpacodiff),axis=2)

newnames=names_aco.copy()
for i in range(len(newnames)):
    newnames[i]=newnames[i]+'Diff'
Ns[0]=np.concatenate((Ns[0],newnames),axis=0)


#tmpacodiff=np.diff(tmpacodiff, n=1, axis=0)

# add crude temporal differential of audio features after downsampling
#Xs{1} = cat(3,Xs{1},Xs{1}-[zeros(size(Xs{1}(1,:,:,:)));Xs{1}(1:end-1,:,:,:)]);  %referencing also 4th dimension for Paul's data
#Ns{1} = cat(1,Ns{1},celfun(@(x)[x,'Diff'],Ns{1})); %add name of differential of audio features
    
    












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



