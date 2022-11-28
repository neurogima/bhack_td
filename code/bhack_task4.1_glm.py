# Do the imports

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





## brainhack task 04 Python translated: GLMs and CV stats

if not do_brainhack or brainhack_task == 4 #do this if not in brainhack mode, or if brainhack task == 3

    if do_brainhack #load the required if in brainhack mode
        print('loading intermediate variables for task 3')
        print(fnout_tasks[2])
        out_dist = io.loadmat(fnout_tasks[2])
        print('done')
    
    ## fit the distance-based GLMs, and generalize to test data
    Ds = out_dist['Ds']
    Y = Ds[0, -1] #iEEG data, of size [ntimepairs nsensors ntimechunks]
    Y = np.transpose( np.expand_dims(Y, axis=1), (0,1,3,2)) #iEEG data, of size [ntimepairs 1 ntimechunks nsensors]
    #we need this shape for the iEEG data matrix to
    #facilitate the GLM modelling, below
    
    # useful inline functions taken from BLG_GLM_ND used to recompute
    # cross-validated predictions of the GLM models
    
    ## define demean function
    
#######################################################
    fdemean=@(x) bsxfun(@minus,x,mean(x)); %demean
    fregrdem=@(x) cat(2,ones(size(x(:,1,:))),fdemean(x)); %demean and add intercept
    
    ## we already compute the total sum of squares of each temporal chunk for
    ## each sensor. It's used to compute the RSQ_cv, below
    SSTtest=sum(bsxfun(@minus,Y,mean(Y)).^2);
    
    
    Stats=cell(0);
    c0=clock; %current time
    disp('fitting GLMs')
    for whichpred=1:2
        
        [BetaTrain1,~,~,~]=BLG_GLM_ND(Ds{whichpred},Y,0,0); %fit one GLM for each temporal chunk, and for each sensor
        %output is matrix Beta of GLM coefficients, of size:
        % [npredictors+1 1 ntemporalchunks nsensors];
        % note that BetaTrain1(1,:,:,:), is the GLM beta for the intercept
        % term, i.e., BetaTrain1(2,:,:,:) is the GLM beta for the first model predictor in Preds{1};
        
        %%% for the moment, we implement a flavour of
        %%% leave-one-temporal-chunk-out cross-validation scheme. In
        %%% particular, we consider as training GLM coefficients for temporal
        %%% chunk 1 the average GLM coefficient of temporal chunks 2:n_chunks
        %%% The particular CV scheme is likely to be revised in future
        %%% iterations. This will however do not require drastic changes to the
        %%% pipeline.
        BetaTrain2=zeros(size(BetaTrain1)); %here we create the matrix of GLM coefficients averaged across training temporal chunks
        
        for ichunk=1:n_chunks %%% loop through temporal chunks
            idx_train_chunks=setxor(ichunk,1:n_chunks); %indices to temporal chunks across which we will average
            %%% the training GLM coefficients. Excludes the current temporal chunk
            BetaTrain2(:,:,ichunk,:)=mean(BetaTrain1(:,:,idx_train_chunks,:),3); %average betas across training temporal chunks
        end
        
        
        
        PredTest=mtimesx(fregrdem(Ds{whichpred}),BetaTrain2); %test-set prediction based on training-set GLM betas
        
        SSEtest=sum(bsxfun(@minus,Y,PredTest).^2,1); %sum of squared errors for test-set prediction
        
        RSQ_cv=1-SSEtest./SSTtest; %cross-validated R squared
        
        r_cv = BLGmx_corr2(Y,PredTest); %cross-validated correlation
        
        Stats{whichpred,1}=RSQ_cv; %add to output cell
        Stats{whichpred,2}=r_cv; %add to output cell
        
        et=etime(clock,c0); %time elapsed since beginning of this four loop
        disp(['---All done for model: ',num2str(whichpred),' elapsed time: ',num2str(et)])
    end
    
    
    tmpstats=celfun(@(x) (prctile(x,50,3)),Stats);
    tmp=cell2mat(tmpstats);
    sensornames={'sensitive to acoustics' 'sensitive to semantics' 'unresponsive'};
    statnames={'RSQ_cv' 'r_cv'};
    modelnames={'aco model' 'sem model'};
    for i=1:length(sensornames)
        disp(['Sensor: ',sensornames{i}])
        
        t=array2table(double(tmp(:,:,i)),'VariableNames',statnames,'RowNames',modelnames);
        disp(t)
        
    end
    
end

if save_intermediate && do_brainhack==0 %save intermediate variables, if requested and not in brainhack mode
    disp('saving intermediate file: ')
    disp(fnout_tasks{4})
    save(fnout_tasks{4},'Stats')
    disp('done')
end


