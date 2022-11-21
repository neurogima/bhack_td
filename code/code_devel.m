
%% start devel

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


%%
%%% set critical analysis parameters
fs = 125; %target sampling frequency for all signals

max_time_s = 200;%maximum time of timecourse to consider for analyses (seconds)

n_chunks = 20; %divide the signal into nchunks "independent" time-courses.
              %GLMs will be cross-validated across these nchunks. For
              %Paul each trial is a chunk

p_chunk_out = 0.05; %proportion of chunk time points that will be removed 
              %from the beginning (50%)/end(50%) of each chunk. Discarding some of the
              %timepoints connecting the different chunks increases the
              %independence of chunks, and betters the generalization.
              %Also potentially useful if chunks are trials and onset/offset
              %effects need to be discarded from analysis
              
ms_lags = 0:8:100; %feature-to-brain-lags (ms; 200 ms = the brain represents
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

ns_time_folding = floor(ms_time_folding/1000*fs); %folding window in n samples



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

clear tmpout tmp time_sem time_eeg time_aco tmpfs tmpdown stmp stmp2 names_sem names_aco names_eeg i j ilag fss fs_aco fs_sem fs_eeg dat_fns X_aco X_sem X_eeg Xs tmpisnan tmpnanmean

if do_zscore %do z-scoring of distances on a chunk-by-chunk basis, if requested.
             %this is recommended for cross-validation at the moment.
    Ds=celfun(@zscore,Ds);
end














