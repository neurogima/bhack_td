%% set paths and analysis parameters: run for all tasks
clear
clc

rootmain='D:/!!Projects/bhack_td/';
cd(rootmain)
run([rootmain,'code/project_init.m'])


do_brainhack=1; %brainhack mode;

brainhack_task = 4; %which brainhack task? 
   %determines which intermediate variables will be loaded

save_intermediate=0; %save intermediate variables

%%% filenames for intermediate variables
fnout_pars=[rootmain,'results/bhack_analysis_pars.mat'];
fnout_tasks={'bhack_task_01_output.mat',...
             'bhack_task_02_output_temporalfolding.mat',...
             'bhack_task_03_output_distances.mat',... %this file is not stored on github, see below
             'bhack_task_04_output_GLMs.mat',...
             'bhack_task_05_output_plots.pdf'};

%add path to out filenames for different tasks        
fnout_tasks=celfun(@(x) [rootmain,'results/',x],fnout_tasks);

if ~exist(fnout_tasks{3},'file')
    disp('bhack_task_03_output_distances is not available in repo')
    disp('download from : ')
    disp('https://drive.google.com/file/d/1b5vRLCH3fnvt4EtXwUoT2JmhDZg9h6Cl/view?usp=sharing')
    disp('after downloading place in results/ folder')
end

%this file is not stored on github (390 MB)
% 'bhack_task_03_output_distances.mat',... 
% download from:
% after downloading place in results folder

%% set analysis parameters: run for all tasks

if ~do_brainhack %if not brainhack mode declare parameters
    
    %%% set critical analysis parameters
    %%% caution: high fs, high max_time_s, low n_chunks and low ms_time_folding
    %%% and many different lags will all result in computationally intensive
    %%% analyses that will likely saturate your RAM and potentially return an
    %%% error
    
    fs = 125/4; %target sampling frequency for all signals
    
    max_time_s = 500;%maximum time of timecourse to consider for analyses (seconds)
    
    n_chunks = 10; %divide the signal into nchunks "independent" time-courses.
    %GLMs will be cross-validated across these nchunks. For
    %Paul each trial is a chunk
    
    p_chunk_out = 0.05; %proportion of chunk time points that will be removed
    %from the beginning (50%)/end(50%) of each chunk. Discarding some of the
    %timepoints connecting the different chunks increases the
    %independence of chunks, and betters the generalization.
    %Also potentially useful if chunks are trials and onset/offset
    %effects need to be discarded from analysis
    
    ms_lags = 100:50:500; %feature-to-brain-lags (ms; 200 ms = the brain represents
    %the feature 200 ms after the waveform reaches the tympanic
    %membrane)
    
    ms_time_folding= 100; %this is the folding trick for decreasing the size of the temporal
    %distance matrix. Instead of computing the distance between
    %timepoint A and timepoint B, we compute the distance between
    %timecourse A and timecourse B, of duration ms_time_fold.
    
    do_zscore = 1; %standardize distances (independently for each chunk)
    %before computing any GLM
    
    
    
    %%% derive analysis parameters in number of samples
    ns_max = max_time_s*fs; %maximum number of samples of the considered timecourse.
    
    ns_chunk = floor(ns_max/n_chunks); %number of samples in each chunk
    
    ns_max = n_chunks*ns_chunk; %recompute max samples, after accounting for rounding
    
    ns_chunk_out = floor(ns_chunk*p_chunk_out); %how many samples of the chunk will be removed?
    %50% from beginning, 50% from end, so
    %ns_chunk_out must be even
    
    ns_chunk_out = floor(ns_chunk_out/2)*2; %make ns_chunk_out even
    
    ns_lags = unique(round(ms_lags./1000.*fs)); %feature-to-brain lags in n samples
    
    n_lags = length(ns_lags); %number of considered lags
    
    ms_lags = ns_lags./fs*1000;
    
    ns_time_folding = floor(ms_time_folding/1000*fs); %folding window in n samples
    
    
    %store analysis parameters
    analysis_pars=struct([]);
    analysis_pars(1).fs = fs;
    analysis_pars.max_time_s = max_time_s;
    analysis_pars.n_chunks = n_chunks;
    analysis_pars.p_chunk_out = p_chunk_out;
    analysis_pars.ms_lags = ms_lags;
    analysis_pars.ms_time_folding= ms_time_folding;
    analysis_pars.do_zscore = do_zscore;
    analysis_pars.ns_max = ns_max;
    analysis_pars.ns_chunk = ns_chunk;
    analysis_pars.ns_chunk_out = ns_chunk_out;
    analysis_pars.ns_lags = ns_lags;
    analysis_pars.n_lags = n_lags;
    analysis_pars.ns_time_folding = ns_time_folding;
    analysis_pars_cell={fs max_time_s n_chunks p_chunk_out ms_lags ms_time_folding,...
        do_zscore ns_max ns_chunk ns_chunk_out ns_lags n_lags ns_time_folding};
    analysis_par_names={'fs' 'max_time_s' 'n_chunks' 'p_chunk_out' 'ms_lags' 'ms_time_folding',...
        'do_zscore' 'ns_max' 'ns_chunk' 'ns_chunk_out' 'ns_lags' 'n_lags' 'ns_time_folding'};
end

if save_intermediate && do_brainhack==0 %save intermediate variables, if requested and not in brainhack mode
    disp('saving intermediate file: ')
    disp(fnout_pars)
    save(fnout_pars,'analysis_pars','analysis_pars_cell','analysis_par_names')
    disp('done')
end




%% brainhack task 01: preliminaries, data massaging

if ~do_brainhack || brainhack_task==1 %do this if not in brainhack mode, or if brainhack task == 1
    
    if do_brainhack %if brainhack mode load analysis parameters
        load(fnout_pars) %#ok<UNRCH>
        for i=1:length(analysis_pars_cell) %allocate to workspace the variables
            tmpstr=[analysis_par_names{i},'=analysis_pars_cell{i};'];
            disp(tmpstr)
            eval(tmpstr)
        end
    end
    
    %%% load the data
    fn_ieeg=[rootmain,'dat/brain.mat']; %ieeg data
    fn_aco=[rootmain,'dat/aco.mat'];    %all acoustic models
    fn_sem=[rootmain,'dat/sem.mat'];   %semantic models
    dat_fns={fn_ieeg fn_aco fn_sem};
    dat_nams={'eeg' 'aco' 'sem'};
    for i=1:length(dat_fns) %load data.
        %%%% The content of each dataset is four variables:
        %%%% X = the signal of size [ntimepoints nembeddingdimensions n_sensors(iEEG)_OR_nmodels(aco/sem)]
        %%%% fs = the sampling frequency in Hz
        %%%% names = the names for the sensors/models (third dimension of X)
        %%%% time = time vector in ms
        tmp=load(dat_fns{i});
        disp('---------------')
        disp(dat_fns{i})
        vn=fieldnames(tmp);
        for j=1:length(vn) %%%% loop through the variables contained in this .mat dataset
            varnam=vn{j};
            varnew=[varnam,'_',dat_nams{i}]; %new name for the variable allocated to workspace - appends datasat name (eeg/aco/sem)
            eval([varnew,'=tmp.',varnam,';'])
            if strcmp(varnam,'X')
                eval([varnew,'=double(tmp.',varnam,');']) %assign to workspace
            end
            eval(['tmps=size(',varnew,');']) %verbose feedback
            disp([varnew,', var of size: ',num2str(tmps)])
        end
    end
    
    Xs = {X_aco X_sem X_eeg}; %datasets: this could be a python list?
    
    Xs{1}(isinf(Xs{1}))=0; %replace inf values in acoustics with zeroes;
    
    %%% replace time samples without embedding with average embedding of other
    %%% time samples. This is a temporary hack.
    tmpisnan=sum(isnan(X_sem),2)>0;
    for i=1:2 %give average semantic embedding for time samples without average embedding.
        %this should be done independently for each chunk
        tmpnanmean=nanmean(X_sem(:,:,i),1);
        tmpisnan=sum(isnan(X_sem(:,:,i)),2)>0;
        Xs{2}(tmpisnan,:,i)=ones(sum(tmpisnan),1)*tmpnanmean;
    end
    
    Xs = celfun(@single,Xs); %halve precision to halve memory requirements and speed up computations
    Ts = {time_aco time_sem time_eeg}; %time vectors in ms
    Ts = celfun(@round,Ts); %round to the nearest millisecond, eliminates very small timing vector gremlins
    Ns = {names_aco names_sem names_eeg}; %model/sensor names (third dimension of Xs)
    fss = {fs_aco fs_sem fs_eeg}; %sampling frequencies
    
    %%% let's do a very crude downsampling of data to the new sampling
    %%% frequency
    for i = 1:length(Xs)
        tmpfs = fss{i}; %the sampling frequency of this signal in Xs{i}
        tmpdown = tmpfs/fs; %down-sampling factor: take one sample every tmpdown samples
        if rem(tmpdown,1)~=0 || tmpdown<1 %if not integer or if target fs > tmpfs (upsampling)
            error('the ratio of the old fs to the new fs must be an integer; target fs must be < old fs')
        end
        
        %%% apply the downsampling;
        Xs{i} = Xs{i}(1:tmpdown:end,:,:,:); %referencing also 4th dimension for Paul's data
        Ts{i} = Ts{i}(1:tmpdown:end);
    end
    % clear fss tmpfs tmpdown X_aco X_sem X_eeg  %clear workspace of junk
    
    %%% add crude temporal differential of audio features after downsampling
    Xs{1} = cat(3,Xs{1},Xs{1}-[zeros(size(Xs{1}(1,:,:,:)));Xs{1}(1:end-1,:,:,:)]);  %referencing also 4th dimension for Paul's data
    Ns{1} = cat(1,Ns{1},celfun(@(x)[x,'Diff'],Ns{1})); %add name of differential of audio features
    
    
    %%% OK, let's create temporal chunks. These will take the role of trials
    %%% across which the GLMs will be generalized. Temporal chunks are
    %%% concatenated along the fourth dimension of the Xs matrices. If Xs
    %%% already include different trials, concatenated along the fourth
    %%% dimension, this part should be commented out.
    Xs=celfun(@(x) x(1:min([length(x) ns_max]),:,:,:),Xs); %takes ns_max temporal samples of entire timecourse
    Ts=celfun(@(x) x(1:min([length(x) ns_max])),Ts);%takes ns_max temporal samples of entire time vector
    if ~isequal(Ts{1},Ts{2}) || ~isequal(Ts{1},Ts{3}) %sanity check
        error('check time vectors: they should be now equal for all Xs datasets')
    end
    
    
    %%% OK, let's first cut the timecourses into chunks. Chunks take the role
    %%% of trials, and are concatenated in the 4th dimension of Xs matrices
    chunk_idx = reshape(1:length(Ts{1}),[ns_chunk n_chunks]); %indices to first dimension
    % of Xs matrices, reshaped so that each chunk is in a different column
    
    for i=1:length(Xs)
        tmpout=[];
        for j=1:n_chunks
            tmpout=cat(4,tmpout,Xs{i}(chunk_idx(:,j),:,:)); %concatenate temporal chunks in 4th dimension
        end
        Xs{i}=tmpout; %replace old Xs matrix without temporal chunking
    end
    
    %%%let's chop out ns_chunk_out/2 from beginning and end. This aids the
    %%%independence of data in different trials.
    Xs = celfun(@(x) x(ns_chunk_out/2+1:end-ns_chunk_out/2,:,:,:),Xs);
    
    %%%let's remove max(ns_lags) from end of eeg (Xs(3)) to account for the
    %%%maximum feature-to-brain lag
    Xs{3}=Xs{3}(max(ns_lags):end,:,:,:);
    
% else %if in brainhack mode
%     if brainhack_task>1
%         load(fnout{1})
%     end
end

if save_intermediate && do_brainhack==0 %save intermediate variables, if requested and not in brainhack mode
    disp('saving intermediate file: ')
    disp(fnout_tasks{1})
    save(fnout_tasks{1},'Xs','Ts','Ns','fss')
    disp('done')
end



%% brainhack task 02: temporal folding of time-courses

if ~do_brainhack || brainhack_task==2 %do this if not in brainhack mode, or if brainhack task == 2
    
    if do_brainhack %load the required if in brainhack mode
        load(fnout_pars) %#ok<UNRCH> %load analysis parameters
        for i=1:length(analysis_pars_cell) %allocate to workspace the variables
            tmpstr=[analysis_par_names{i},'=analysis_pars_cell{i};'];
            disp(tmpstr)
            eval(tmpstr)
        end
        disp('loading intermediate variables for task 2')
        disp(fnout_tasks{1})
        load(fnout_tasks{1})
        disp('done')
    end
    
    %%% let's do the temporal folding, so that we compute the distance not
    %%% between time-points, but between small segments of the temporal chunk.
    if ns_time_folding > 0 %if temporal folding has been requested do it
        
        n_time_foldings=floor(size(Xs{3},1)/ns_time_folding); %number of "foldings"
        
        %%%indices to the foldings of the temporal signal
        foldings_idx=reshape(1:ns_time_folding*n_time_foldings,[ns_time_folding n_time_foldings])';
        
        
        %%% VERY CUMBERSOME APPROACH for folding iEEG, but it shows the
        %%% mechanics clearly
        %     tmpout=[];
        %     for i=1:n_time_foldings
        %         tmpout=cat(1,tmpout,permute(Xs{3}(foldings_idx(i,:),:,:,:),[2 1 3 4])); %the ns_time_folding samples are concatenated along the second dimension
        %     end
        
        %%% less cumbersome way to fold iEEG
        disp('folding iEEG in time')
        
        tic %tic toc gives the time between occurred since tic at toc
        
        tmp=Xs{3}; %this is the iEEG data
        
        %concatenate the fold indices in one single column, then add a
        %singleton dimension after the first for the subsequent reshaping
        tmp=permute(tmp(1:max(foldings_idx(:)),:,:,:),[1 5 2 3 4]);
        
        stmp=size(tmp); %current size of matrix
        
        %%%reshape so that the first two dimensions contain the "temporal fold"
        %%%in first, and the number of temporal folds in second dimension
        tmp=reshape(tmp,[ns_time_folding n_time_foldings stmp(3:end)]);
        
        %%% permute matrix shape so that first dimension is second and second
        %%% is first. Also remove the leftover singleton dimension 3 that we
        %%% had introduced before
        tmp=permute(tmp,[2 1 4 5 3]);
        
        Xs{3}=tmp; %put back temporally folded iEEG data inside the Xs data cell
        disp('iEEG folded')
        toc
        
        %%% let's then fold the acoustic and semantic time-courses
        %%% and also add the lags (to be concatenated as additional models in
        %%% third dimension).
        
        %     %%% VERY CUMBERSOME APPROACH for folding and lagging acoustics
        %     %%% and semantic models, but shows the mechanics
        %     disp('folding models in time, and adding lags')
        %     for acosem=1:2
        %         tmpout2=[];
        %         for ilag=1:n_lags
        %             tmpout1=[];
        %             for i=1:n_time_foldings
        %                 tmpout1=cat(1,tmpout1,permute(Xs{acosem}(foldings_idx(i,:)+ns_lags(ilag),:,:,:),[2 1 3 4]));
        %                 %%%%the ns_time_folding samples are concatenated along the second dimension
        %                 disp(num2str([i acosem ilag]))
        %             end
        %             tmpout2=cat(3,tmpout2,tmpout1);
        %             disp(num2str([acosem ilag]))
        %         end
        %         Xs{acosem}=tmpout2;
        %     end
        
        %%% less cumbersome way to apply temporal folding, and feature-to-brain
        %%% lagging to the acoustics and semantics signals
        disp('folding acoustics and semantics in time')
        tic
        for acosem=1:2
            tmpout=[];
            for ilag=1:n_lags
                ns_thislag=ns_lags(ilag);
                tmp=Xs{acosem};
                
                %%% select the data considering the temporal folding indices,
                %%% but increase the indices by ns_thislag samples (samples for
                %%% this particular lagged version of the model) so as to
                %%% create the lagged models. The rest of the
                %%% reshaping/permuting work is the same as for the iEEG data,
                %%% above.
                tmp=permute(tmp([1:max(foldings_idx(:))]+ns_thislag,:,:,:),[1 5 2 3 4]); %#ok<NBRAK>
                stmp=size(tmp);
                tmp=reshape(tmp,[ns_time_folding n_time_foldings stmp(3:end)]);
                tmp=permute(tmp,[2 1 3 4 5]);
                stmp=size(tmp);
                tmp=reshape(tmp,[stmp(1) stmp(2)*stmp(3) stmp(4:5)]);
                if ilag==1 %if this is the first lag, preallocate matrix to speed up
                    stmp=size(tmp);
                    stmp2=stmp;
                    stmp2(3)=stmp2(3)*n_lags;
                    tmpout=zeros(stmp2,'single'); %preallocate output matrix to speed up
                end
                
                %%% concatenate lagged versions of the model along the third
                %%% dimension of the data matrix (tmpout).
                tmpout(:,:,[1:size(tmp,3)]+(ilag-1)*size(tmp,3),:)=tmp; %#ok<NBRAK,SAGROW>
                %disp(num2str([acosem ilag]))
            end
            %         tmp2=tmpout;
            Xs{acosem}=tmpout;
        end
        disp('acoustics and semantics folded')
        toc
        
        
    end
end

if save_intermediate && do_brainhack==0 %save intermediate variables, if requested and not in brainhack mode
    disp('saving intermediate file: ')
    disp(fnout_tasks{2})
    save(fnout_tasks{2},'Xs','Ts','Ns','fss')
    disp('done')
end


%% brainhack task 03: distance computation

if ~do_brainhack || brainhack_task==3 %do this if not in brainhack mode, or if brainhack task == 3
    
    if do_brainhack %load the required if in brainhack mode
        load(fnout_pars) %#ok<UNRCH> %load analysis parameters
        for i=1:length(analysis_pars_cell) %allocate to workspace the variables
            tmpstr=[analysis_par_names{i},'=analysis_pars_cell{i};'];
            disp(tmpstr)
            eval(tmpstr)
        end
        disp('loading intermediate variables for task 3')
        disp(fnout_tasks{2})
        load(fnout_tasks{2})
        disp('done')
    end
    
    
    %%% let's compute the distances
    disp('Computing distances')
    tic;
    Ds=cell(0);
    Ds{1}=BLG_CosDistND(permute(Xs{1},[2 1 3 4])); %put features in second dimension and time windows in first dimension before computing distance
    Ds{2}=BLG_CosDistND(permute(Xs{2},[2 1 3 4])); %put features in second dimension and time windows in first dimension before computing distance
    Ds{3}=BLG_EucDistND(permute(Xs{3},[2 1 3 4])); %put features in second dimension and time windows in first dimension before computing distance
    disp('Distances computed')
    toc
    
    %%% deal with NaNs and Infs in distance matrix
    for i=1:length(Ds) %this absolutely ignores the root of the problem, but let's ignore inf and nan values in distances, for the moment
        %this is a very ugly temporary hack.
        Ds{i}(isinf(Ds{i}))=0;
        Ds{i}(isnan(Ds{i}))=nanmean(Ds{i}(:));
    end
    
    
    % if do_zscore %do z-scoring of distances on a chunk-by-chunk basis, if requested.
    %     %this is recommended for cross-validation at the moment.
    %     Ds=celfun(@zscore,Ds);
    % end
    
end


if save_intermediate && do_brainhack==0 %save intermediate variables, if requested and not in brainhack mode
    disp('saving intermediate file: ')
    disp(fnout_tasks{3})
    save(fnout_tasks{3},'Ds','Ts','Ns','fss')
    disp('done')
end


%% brainhack task 04: GLMs and CV stats

if ~do_brainhack || brainhack_task==4 %do this if not in brainhack mode, or if brainhack task == 3
    
    if do_brainhack %load the required if in brainhack mode
        load(fnout_pars) %#ok<UNRCH> %load analysis parameters
        for i=1:length(analysis_pars_cell) %allocate to workspace the variables
            tmpstr=[analysis_par_names{i},'=analysis_pars_cell{i};'];
            disp(tmpstr)
            eval(tmpstr)
        end
        disp('loading intermediate variables for task 3')
        disp(fnout_tasks{3})
        load(fnout_tasks{3})
        disp('done')
    end
    
    %%% fit the distance-based GLMs, and generalize to test data
    
    Y=Ds{3}; %iEEG data, of size [ntimepairs nsensors ntimechunks]
    Y=permute(Y,[1 4 3 2]); %iEEG data, of size [ntimepairs 1 ntimechunks nsensors]
    %we need this shape for the iEEG data matrix to
    %facilitate the GLM modelling, below
    
    % useful inline functions taken from BLG_GLM_ND used to recompute
    % cross-validated predictions of the GLM models
    fdemean=@(x) bsxfun(@minus,x,mean(x)); %demean
    fregrdem=@(x) cat(2,ones(size(x(:,1,:))),fdemean(x)); %demean and add intercept
    
    %%% we already compute the total sum of squares of each temporal chunk for
    %%% each sensor. It's used to compute the RSQ_cv, below
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




%% brainhack task 05: Plots

if ~do_brainhack || brainhack_task==5 %do this if not in brainhack mode, or if brainhack task == 3
    
    if do_brainhack %load the required if in brainhack mode
        
        disp('loading intermediate variables for task 4') %#ok<UNRCH>
        disp(fnout_tasks{4})
        load(fnout_tasks{4})
        disp('done')
    end
    
    %%% some plots
    % close all
    figure
    hold
    whichstat=1;
    ylabs={'R^2_C_V' 'r_C_V'};
    cols=[[1 0 0];[0 0 1]];
    tmp=cell2mat(Stats(:,whichstat));
    tmp=permute(tmp,[3 4 1 2]);
    size(tmp);
    oris={'left' 'right'};
    addx=[0 .4]-0.375;
    for i=[1:2]
        distributionPlot(tmp(:,1,i),'histOri',oris{i},'color',cols(i,:),'widthDiv',[2 2],...
            'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
            'xValues',[1]+addx(i));
    end
    l=legend({'acoustic models' 'semantic models'},'location','northeast','autoupdate','off');
    for i=[2:-1:1]
        distributionPlot(tmp(:,:,i),'histOri',oris{i},'color',cols(i,:),'widthDiv',[2 2],...
            'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
            'xValues',[1:3]+addx(i));
    end
    set(gca,'xtick',1:3,'view',[90 90],'xticklabel',{'acoustic' 'semantic' 'non-responsive'})
    ylabel(ylabs{whichstat})
    xlabel('Sensor type')
    axis tight
    title('Distribution across CV folds')
end


if save_intermediate && do_brainhack==0 %save intermediate variables, if requested and not in brainhack mode
    disp('saving intermediate file: ')
    disp(fnout_tasks{5})
    saveas(gcf,fnout_tasks{5})
    disp('done')
end
















