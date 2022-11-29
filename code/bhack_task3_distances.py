# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 10:40:57 2022

@author: marieplgt
"""
	
# In[set environment]:
    
import scipy as sc
from scipy import io
import numpy as np
import pandas as pd
import h5py
import seaborn
import matplotlib
import nibabel
import os

# In[set paths etc.]:
workdir='D:\Coding\GitHub_Repo/bhack_td/'
out_dir=workdir+'results/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

fname=out_dir+'bhack_task_02_output_temporalfolding.mat'
tmp=io.loadmat(fname)

Xs = tmp['Xs'][0]#halve precision to halve memory requirements and speed up computations
Ts = tmp['Ts'][0] #time vectors in ms
Ns=tmp['Ns'][0] #model/sensor names (third dimension of Xs)
fss = tmp['fss'][0] #sampling frequencies

# In[Translation]:

ac_data=np.asarray(Xs[0])
sem_data=Xs[1]
eeg_data=Xs[2]


stmp=ac_data.shape
ntimepairs=stmp[0]*(stmp[0]-1)/2
snew=[int(ntimepairs),stmp[2],stmp[3]]
for i in range(sizetmp[2]):
    for j in range(sizetmp[3]):
        Ds_ac[:,i,j] = sc.spatial.distance.pdist(ac_data[:,:,i,j], metric='cosine')
        print([i,j])
        Xs[0]=Ds_ac

stmp2=sem_data.shape
snew2=[int(ntimepairs),stmp2[2],stmp2[3]]
for i in range(stmp2[2]):
    for j in range(stmp2[3]):
        Ds_sem[:,i,j] = sc.spatial.distance.pdist(sem_data[:,:,i,j], metric='cosine')
        print([i,j])
        Xs[1]=Ds_sem


stmp3=eeg_data.shape
snew3=[int(ntimepairs),stmp3[2],stmp3[3]]
for i in range(stmp3[2]):
    for j in range(stmp3[3]):
        Ds_eeg[:,i,j] = sc.spatial.distance.pdist(eeg_data[:,:,i,j], metric='euclidean')
        print([i,j])
        Xs[2]=Ds_eeg

print(Ds_ac.shape)
Xs[0].shape
print(Ds_sem.shape)
Xs[1].shape
print(Ds_eeg.shape)
Xs[2].shape

















    
    