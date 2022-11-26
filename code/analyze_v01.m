%% set paths, load data. Modify rootmain as required to set local paths

rootmain='D:/!!Projects/bhack_td/';
% do_ft=1;
cd(rootmain)
run([rootmain,'code/project_init.m'])

fn_ieeg=[rootmain,'dat/brain.mat']; %ieeg data
fn_aco=[rootmain,'dat/aco.mat'];    %all acoustic models
fn_sem=[rootmain,'dat/sem.mat'];   %semantic models



dat_fns={fn_ieeg fn_aco fn_sem};
dat_nams={'eeg' 'aco' 'sem'};
for i=1:length(dat_fns)
    tmp=load(dat_fns{i});
    disp('---------------')
    disp(dat_fns{i})
    vn=fieldnames(tmp);
    for j=1:length(vn)
        varnam=vn{j};
        varnew=[varnam,'_',dat_nams{i}];
        eval([varnew,'=tmp.',varnam,';'])
        if strcmp(varnam,'X')
            eval([varnew,'=double(tmp.',varnam,');'])
        end
        eval(['tmps=size(',varnew,');'])
        disp([varnew,', var of size: ',num2str(tmps)])
    end
end

clear varnew varnam vn tmp i j RESTOREDEFAULTPATH_EXECUTED tmps %clear workspace of junk
% %%% loaded variables:
% rootmain,'/dat/brain.mat'
% X_eeg, var of size: 288350       1       3 % sensors in third dimension
% fs_eeg, var of size: 1  1 %fs of iEEG data = 500 Hz
% names_eeg, var of size: 3  1 %name of iEEG sensors
% time_eeg, var of size: 288350       1 %time vector of iEEG data in ms
% ---------------
% rootmain,'/dat/aco.mat'
% X_aco, var of size: 576499       1       5 %acoustic models in third dimension
% fs_aco, var of size: 1  1 %fs of acoustic models = 1000 Hz
% names_aco, var of size: 5  1 %name of acoustic models
% time_aco, var of size: 576499       1 %time vector of acoustic models in ms
% ---------------
% rootmain,'/dat/sem.mat'
% X_sem, var of size: 300001     512       2 %semantic models in third dimension
% fs_sem, var of size: 1  1 %fs of semantic model = 500 Hz
% names_sem, var of size: 2  1 %name of semantic models
% time_sem, var of size: 300001       1 %time vector of semantic models in ms


%% set analysis parameters, and compute the distances
%%%
%%%
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
analysis_pars_cell={fs max_time_s n_chunks p_chunk_out ms_lags ms_time_folding,...
    do_zscore ns_max ns_chunk ns_chunk_out ns_lags n_lags}; 
analysis_par_names={'fs' 'max_time_s' 'n_chunks' 'p_chunk_out' 'ms_lags' 'ms_time_folding',...
    'do_zscore' 'ns_max' 'ns_chunk' 'ns_chunk_out' 'ns_lags' 'n_lags'};




Xs = {X_aco X_sem X_eeg}; %datasets

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

%%% clear junk from workspace
% clear tmpout tmp time_sem time_eeg time_aco tmpfs tmpdown stmp stmp2 names_sem names_aco names_eeg i j ilag fss fs_aco fs_sem fs_eeg dat_fns X_aco X_sem X_eeg Xs tmpisnan tmpnanmean

if do_zscore %do z-scoring of distances on a chunk-by-chunk basis, if requested.
    %this is recommended for cross-validation at the moment.
    Ds=celfun(@zscore,Ds);
end

%% fit the distance-based GLMs, and generalize to test data

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

%% some plots
% close all
figure
hold
whichstat=1;
ylabs={'R^2_C_V' 'r_C_V'};
cols=[[1 0 0];[0 0 1]];
tmp=cell2mat(Stats(:,whichstat));
tmp=permute(tmp,[3 4 1 2]);
size(tmp)
oris={'left' 'right'}
addx=[0 .4]-0.375
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



%% OK, let's do the variance partitioning to reveal the unique contribution of acoustics
%% and semantics, as well as the common acoustics/semantics contribution to
%% the prediction of the temporal iEEG distance

%%% for the computation of the variance partitions according to the
%%% commonality analysis equations, we need to fit three models:
%%% 1. acoustics only
%%% 2. semantics only
%%% 3. acoustics+semantics
%%% Above, we have analysed models 1. and 2. We now need to compute the CV
%%% stats for model 3.



%% fit the distance-based GLM for aco+sem model, and generalize to test data

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

X=cat(2,Ds{1},Ds{2}); %concatenate acoustic and semantic predictors along
     % the second dimension of the predictors matrix X




Stats_model3=cell(0);
c0=clock; %current time

[BetaTrain1,~,~,~]=BLG_GLM_ND(X,Y,0,0); %fit one GLM for each temporal chunk, and for each sensor
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



PredTest=mtimesx(fregrdem(X),BetaTrain2); %test-set prediction based on training-set GLM betas

SSEtest=sum(bsxfun(@minus,Y,PredTest).^2,1); %sum of squared errors for test-set prediction

RSQ_cv=1-SSEtest./SSTtest; %cross-validated R squared

r_cv = BLGmx_corr2(Y,PredTest); %cross-validated correlation

Stats_model3{1,1}=RSQ_cv; %add to output cell
Stats_model3{1,2}=r_cv; %add to output cell

et=etime(clock,c0); %time elapsed since beginning of this four loop
disp(['---All done for aco+sem model; elapsed time: ',num2str(et)])


%% do variance partitioning
Stats3=cat(1,Stats(:,1),Stats_model3(:,1)); %aggregate RSQ_CV for models 1/2/3
Stats3=celfun(@squeeze,Stats3); %squeeze matrices, i.e. remove leading singleton dimensions of the matrix (e.g., first two dimensions of matrix tmp of size [1 1 3 4])
Stats3_m=cell2mat(permute(Stats3,[3 4 1 2])); %convert cell to matrix to facilitate coding
Stats3_m=permute(Stats3_m,[2 3 1]); %put temporal chunks last dimension
size(Stats3_m)
%%% Stats3_m is of size [nsensors=[aco sem nonresp],
%%% n_models=[models1,2,3], n_temporal_chunks];

%%% do variance partitioning
VarPart3_m=zeros(size(Stats3_m));

%unique(aco) = RSQ_cv(aco+sem) - RSQcv(sem)
VarPart3_m(:,1,:)=Stats3_m(:,3,:)-Stats3_m(:,2,:); 

%unique(sem) = RSQ_cv(aco+sem) - RSQcv(aco)
VarPart3_m(:,2,:)=Stats3_m(:,3,:)-Stats3_m(:,1,:);

%common(aco,sem) = RSQ_cv(aco) + RSQ_cv(sem) - RSQ_cv(aco+sem)
VarPart3_m(:,3,:)=Stats3_m(:,1,:)+Stats3_m(:,2,:)-Stats3_m(:,3,:);

% VarPart=cell(size(Stats3));

VarPart3=cat(1,Stats(:,:),Stats_model3(:,:)); %aggregate RSQ_CV for models 1/2/3


%% some plots
% close all
figure
hold
whichstat=1;
ylabs={'R^2_C_V' 'r_C_V'};
cols=[[1 0 0];[0 0 1];[0 1 0]];
tmp=permute(VarPart3_m,[3 1 2]);
size(tmp)
oris={'left' 'right' 'left'}
addx=[-0.4 0 -0.4]-0.375
for i=[1:3]
    distributionPlot(tmp(:,1,i),'histOri',oris{i},'color',cols(i,:),'widthDiv',[2 2],...
        'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
        'xValues',[1]+addx(i));
end

l=legend({'unique acoustics' 'unique semantics' 'common acoustics semantics'},'location','northeast','autoupdate','off');
for i=[3:-1:1]
    distributionPlot(tmp(:,:,i),'histOri',oris{i},'color',cols(i,:),'widthDiv',[2 2],...
        'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
        'xValues',[1:3]+addx(i));
end

set(gca,'xtick',1:3,'view',[90 90],'xticklabel',{'acoustic' 'semantic' 'non-responsive'})
ylabel(ylabs{whichstat})
xlabel('Sensor type')
axis tight
title('Variance partitioning: Distribution across CV folds')


%% do models for each different feature-to-brain lag

Ds_bylag=cell(0); %lag-specific distances
for i=1:2 %loop through acoustic and semantic models
    n_features=size(Ds{i},2)/n_lags;
    for j=1:n_lags %loop through lags
        idx=((j-1)*n_features+1):j*n_features;
        Ds_bylag{j,i}=Ds{i}(:,idx,:);
    end
end




%% fit the distance-based GLMs, and generalize to test data

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


Stats_bylag=cell(0);
c0=clock; %current time
disp('fitting GLMs')
for whichpred=1:2
    for whichlag=1:n_lags
    
    [BetaTrain1,~,~,~]=BLG_GLM_ND(Ds_bylag{whichlag,whichpred},Y,0,0); %fit one GLM for each temporal chunk, and for each sensor
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
    
    
    
    PredTest=mtimesx(fregrdem(Ds_bylag{whichlag,whichpred}),BetaTrain2); %test-set prediction based on training-set GLM betas
    
    SSEtest=sum(bsxfun(@minus,Y,PredTest).^2,1); %sum of squared errors for test-set prediction
    
    RSQ_cv=1-SSEtest./SSTtest; %cross-validated R squared
    
    r_cv = BLGmx_corr2(Y,PredTest); %cross-validated correlation
    
    Stats_bylag{whichlag,whichpred,1}=RSQ_cv; %add to output cell
    Stats_bylag{whichlag,whichpred,2}=r_cv; %add to output cell
    
    et=etime(clock,c0); %time elapsed since beginning of this four loop
    disp(['---All done for model: ',num2str([whichpred whichlag]),' elapsed time: ',num2str(et)])
    end
end

%% some plots
close all
tmprsqcv=cell2mat(Stats_bylag(:,:,1));
tmprsqcv=permute(tmprsqcv,[3 1 2 4]); %n_chunks n_lags n_models n_sensors
size(tmprsqcv)


sensnam={'Sensor: sensitive to acoustic' 'Sensor: sensitive to semantic' 'Sensor: non-responsive'};
figure

whichstat=1;
ylabs={'R^2_C_V' 'r_C_V'};
cols=[[1 0 0];[0 0 1];[0 1 0]];
oris={'left' 'right' 'left'}
addx=[-0.4 0 -0.4]
for whichsens=1:3
    subplot(3,1,whichsens)
    
    hold
    
    if whichsens==1
        for whichmodel=1:2
            tmp=tmprsqcv(:,1,whichmodel,whichsens);
            distributionPlot(tmp,'histOri',oris{whichmodel},'color',cols(whichmodel,:),'widthDiv',[2 2],...
                'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
                'xValues',[1]+addx(whichmodel));
        end
        l=legend({'acoustics model' 'semantics model'},'location','northwest','autoupdate','off');
        
    end
    
    for whichmodel=1:2
        tmp=tmprsqcv(:,:,whichmodel,whichsens);
        distributionPlot(tmp,'histOri',oris{whichmodel},'color',cols(whichmodel,:),'widthDiv',[2 2],...
        'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
        'xValues',[1:n_lags]+addx(whichmodel));
    end
    xt=get(gca,'xtick');
    set(gca,'xtick',1:n_lags,'xticklabel',num2str(ms_lags'),'xticklabelrotation',45)
    
    title(sensnam{whichsens})
    if whichsens==3
        xlabel('Feature-to-brain lag (ms)')
    end
    
    ylabel(ylabs{whichstat})
    plot(get(gca,'xlim'),[0 0],'k--')
    axis tight
end

%%
oris={'left' 'right' 'left'}
addx=[-0.4 0 -0.4]-0.375
for i=[1:2]
    distributionPlot(tmp(:,1,i),'histOri',oris{i},'color',cols(i,:),'widthDiv',[2 2],...
        'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
        'xValues',[1]+addx(i));
end

l=legend({'unique acoustics' 'unique semantics' 'common acoustics semantics'},'location','northeast','autoupdate','off');
for i=[3:-1:1]
    distributionPlot(tmp(:,:,i),'histOri',oris{i},'color',cols(i,:),'widthDiv',[2 2],...
        'addBoxes',0,'showMM',0,'globalNorm',0,'distWidth',0.75,...
        'xValues',[1:3]+addx(i));
end

set(gca,'xtick',1:3,'view',[90 90],'xticklabel',{'acoustic' 'semantic' 'non-responsive'})
ylabel(ylabs{whichstat})
xlabel('Sensor type')
axis tight
title('Variance partitioning: Distribution across CV folds')






