
%% select a couple of sensors for brainhack
rootmain='D:/!!Projects/DFI_iEEG/';
cd(rootmain)
do_ft=0;
run([rootmain,'code/project_init.m'])


Subjs={'subj-01'
    'subj-02'
    'subj-03'
    'subj-04'
    'subj-05'
    'subj-06'
    'subj-07'
    'subj-08'
    'subj-09'
    'subj-10'
    'subj-11'
    'subj-12'
    'subj-13'
    'subj-14'
    'subj-15'
    'subj-16'
    'subj-17'
    'subj-18'};
nsubjs=length(Subjs);



downsampfactor=15; %take one every downsampfactor samples of the 500 Hz timecourse
do_orth=0;
nswin=30;
dodist=1;
doabs=0;

downsampfactor=15; %take one every downsampfactor samples of the 500 Hz timecourse
do_orth=0;
nswin=30;
dodist=1;
lags_ms=[0:20:300];
nlags=length(lags_ms);


downsampfactor=20; %take one every downsampfactor samples of the 500 Hz timecourse
do_orth=1;
nswin=45;
dodist=1;
lags_ms=[0:20:300];
nlags=length(lags_ms);

dur_slice_ms=nswin*1000;
if dodist==0
    adddistfn='nodist';
elseif dodist==1
    adddistfn='dist';
end
addfnam=['_slice_allsamp_',num2str(dur_slice_ms/1000),'s_win_',num2str(downsampfactor),'_downsamp_',num2str(nlags),'_nlags_',adddistfn];
if do_orth==0
    addfnam=['_corrs',addfnam];
elseif do_orth==1
    addfnam=['_pcorrs',addfnam];
end
fnsdatcorr=celfun(@(x) [rootmain,'participants/',x,'/stats/lagged_NLP',addfnam,'.mat'],Subjs);

fnsdatcorr=celfun(@(x) [rootmain,'participants/',x,'/stats/lagged',addfnam,'_aud_res_nlp.mat'],Subjs);


% fnsdat=celfun(@(x) [rootmain,'participants/',x,'/raw.fif'],Subjs);
fnselec=celfun(@(x) [rootmain,'participants/',x,'/electrodes_info.mat'],Subjs);

electrode_info=celfun(@load,fnselec);


template_fn=[rootmain,'group/icbm152_no_cereb_graynew.nii'];
% v=spm_vol(template_fn);
% [Y,XYZ]=spm_read_vols(v);
% %%
[newbb,newidx]=BLG_worldbb_thr(0.01,template_fn);
%    -73  -108   -66
%     73    75    84


targetvoxdim=[1 1 1];
% resize_img_BLG_nearNeighb([rootmain,'group/icbm152_no_cereb_graynew_222vox.nii'], targetvoxdim, newbb, 1);
template_fn2=[rootmain,'group/ricbm152_no_cereb_graynew_111vox.nii'];
% template_fn2=[rootmain,'group/ricbm152_no_cereb_graynew.nii'];
v=spm_vol(template_fn2);
[Y,XYZ]=spm_read_vols(v);

% targetvoxdim=ones(1,3)*gridres*10;
% BB=nan(2,3); %if nans it uses bounding box of original volume
% resize_img_BLG_nearNeighb([outdir,maskfnam2], targetvoxdim, BB, 1);



fnout=celfun(@(x) [rootmain,'participants/',x,'/electrodes_MNIgridmapping111.mat'],Subjs);
for i=1:nsubjs
    if ~exist(fnout{i},'file')
        thismni=electrode_info{i}.electrode_xyz_mni;
        nsens=size(thismni,1);
        out_idx=zeros([nsens,1]);
        for j=1:nsens
            thisXYZ=bsxfun(@minus,XYZ,thismni(j,:)');
            thisXYZ=sum(thisXYZ.^2);
            tmp=find(thisXYZ==min(thisXYZ));
            tmp=tmp(1);
            out_idx(j)=tmp;
            disp(num2str([i j]))
        end
        out_mni=XYZ(:,out_idx)';
        save(fnout{i},'out_idx','out_mni')
    end
    disp(fnout{i})
end

% %%% let's do some stats
dat=celfun(@load,fnout);


if dodist==0
    adddistfn='nodist';
elseif dodist==1
    adddistfn='dist';
end
addstr=['_slice_allsamp_',num2str(dur_slice_ms/1000),'sw_',num2str(downsampfactor),'ds_',num2str(nlags),'_nlags_',adddistfn];
if do_orth==0
    addstr=['_corrs',addstr];
elseif do_orth==1
    addstr=['_pcorrs',addstr];
end

if doabs
    addstr2=[addstr,'Abs'];
else
    addstr2=[addstr,'Raw'];
end




datcorr=celfun(@load,fnsdatcorr);


datout=[];

for i=1:length(datcorr)
    
    tmpc=datcorr{i}.c;
    tmpc=BLG_FisherZ(tmpc);
    tmpc=squeeze(mean(tmpc,2))'; %average across segments
    [nsens,nfeat]=size(tmpc);
    head=[i*ones(nsens,1) [1:nsens]'];
    datout=cat(1,datout,[head tmpc]);
end
close all
tmp=[datout(:,3:4) datout(:,3)-datout(:,4) datout(:,4)-datout(:,3)];
figure,hist(tmp,50);

prc=80;
idxs=[];
nams={'aco' 'sem' 'acosem' 'semaco'};
prcs=prctile(tmp,prc);
x=abs(bsxfun(@minus,tmp,prc));
mx=min(x);
xm=bsxfun(@minus,x,mx);

%%% OK, select sensors 881 for max aco (and contrast), sensor 521 for max sem (and contrast), 
%%% and sensor 2715 for minimum effect

do_ft=1;
run([rootmain,'code/project_init.m'])

idxs=[881 521 2715];
finalsel=[datout(idxs,:) tmp(idxs,3:4)];
disp(num2str(finalsel))

fnsdat=celfun(@(x) [rootmain,'participants/',x,'/raw.fif'],Subjs);

outcell=cell(0);
for i=1:length(idxs)
    isubj=finalsel(i,1);
    isens=finalsel(i,2);
    
    warning off
    disp(['reading data for participant ',num2str([isubj])])
    dat=ft_read_data(fnsdat{isubj})';
    
    head=ft_read_header(fnsdat{isubj});
    clc
    disp([' data for participant read ',num2str([isubj])])
    disp('data read')
    warning on
    dat=single(dat);
    [nt_eeg,~]=size(dat);
    t_eeg=[0:nt_eeg-1]'./head.Fs.*1000;
    %     idxout=find(t_eeg<max(lags_ms));
    %     t_eeg(idxout)=[];
    %     dat(idxout,:)=[];
    [nt_eeg,nsens]=size(dat);
    
    outcell{i,1}=t_eeg;
    outcell{i,2}=double(dat(:,isens));
    outcell{i,3}=head.Fs;
    clc
    disp([' data for participant read ',num2str([isubj])])
    disp('data read')
    
end

% dat_ieeg=
rootmain='D:/!!Projects/DFI_iEEG/';
cd(rootmain)
do_ft=0;
run([rootmain,'code/project_init.m'])



if ~isequal(outcell{1,1},outcell{3,1}) || ~isequal(outcell{1,1},outcell{2,1}) 
    error('check consistency of time dimensions')
end

%%%
% fnout=([rootmain,'group/ieeg_sel_brainhack.mat']);
fnout=['D:/!!Projects/bhack_td/dat/brain.mat'];
time=outcell{i,1};
X=cell2mat(outcell(:,2)');
X=permute(X,[1 3 2]); %sensor is third dimension
fs=outcell{1,3};
names={'iEEG1' 'iEEG2' 'iEEG3'}';
save(fnout,'X','names','time','fs','-v7.3')
disp(fnout)

%% select semantic models data
clear
rootmain='D:/!!Projects/DFI_iEEG/';
cd(rootmain)
do_ft=0;
run([rootmain,'code/project_init.m'])

fnin=[rootmain,'stimuli/NLPmodels_timevarying.mat'];
fnout=['D:/!!Projects/bhack_td/dat/sem.mat'];

tmpdat=load(fnin);
disp('dat loaded')

% save(fnout,'time_ieeg','dat_ieeg','fs_ieeg','-v7.3')
% disp(fnout)
X=cat(3,tmpdat.dat{1,2},tmpdat.dat{2,2});
time=tmpdat.time'*1000;
names=tmpdat.names(1:2);
fs=tmpdat.fs;
save(fnout,'X','names','time','fs','-v7.3')
disp(fnout)

%% select acoustic models data
clear
rootmain='D:/!!Projects/DFI_iEEG/';
cd(rootmain)
do_ft=0;
run([rootmain,'code/project_init.m'])

fnin=[rootmain,'stimuli/AuditoryDimensions_Clean.mat'];
fnout=['D:/!!Projects/bhack_td/dat/aco.mat'];

tmpdat=load(fnin);
disp('dat loaded')

% save(fnout,'time_ieeg','dat_ieeg','fs_ieeg','-v7.3')
% disp(fnout)
X=permute(tmpdat.X,[1 3 2]); %model in third dimension
time=tmpdat.time;
names=tmpdat.names';
fs=tmpdat.fs;
save(fnout,'X','names','time','fs','-v7.3')
disp(fnout)


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
%%% time samples
tmpisnan=sum(isnan(X_sem),2)>0;
for i=1:2 %give average semantic embedding for time samples without average embedding.
         %this should be done independently for each chunk
    tmpnanmean=nanmean(X_sem(:,:,i),1);
    tmpisnan=sum(isnan(X_sem(:,:,i)),2)>0;
    Xs{2}(tmpisnan,:,i)=ones(sum(tmpisnan),1)*tmpnanmean;
end

Xs = celfun(@single,Xs); %halve precision to halve memory requirements and speed up
Ts = {time_aco time_sem time_eeg}; %time vectors
Ts = celfun(@round,Ts); %round to the nearest millisecond, eliminates very small timing vector gremlins
Ns = {names_aco names_sem names_eeg}; %model/sensor names (third dimension of Xs)
fss = {fs_aco fs_sem fs_eeg}; %sampling frequencies

%%% let's do a very crude downsampling of data to the new sampling
%%% frequency
for i = 1:length(Xs)
    tmpfs = fss{i};
    tmpdown = tmpfs/fs;
    if rem(tmpdown,1)~=0 || tmpdown<1 %if not integer or if target fs > tmpfs (upsampling)
        error('the ratio of the old fs to the new fs must be an integer; target fs must be < old fs')
    end
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

chunk_idx = reshape(1:length(Ts{1}),[ns_chunk n_chunks]); %indices to first dimension
               % of Xs matrices, reshaped so that each chunk is in a different column

for i=1:length(Xs)
    tmpout=[];
    for j=1:n_chunks
        tmpout=cat(4,tmpout,Xs{i}(chunk_idx(:,j),:,:));
    end
    Xs{i}=tmpout;
end

%%%let's chop out ns_chunk_out/2 from beginning and end
Xs = celfun(@(x) x(ns_chunk_out/2+1:end-ns_chunk_out/2,:,:,:),Xs);

%%%let's remove max(ns_lags) from end of eeg (Xs(3)) to account for the
%%%maximum feature-to-brain lag
Xs{3}=Xs{3}(max(ns_lags):end,:,:,:);

%%% let's do the folding
if ns_time_folding > 0
    n_time_foldings=floor(size(Xs{3},1)/ns_time_folding); %number of "foldings"
    foldings_idx=reshape(1:ns_time_folding*n_time_foldings,[ns_time_folding n_time_foldings])'; %indices
                     %%% to the foldings of the temporal signal
    
    %%% let's first fold the iEEG signal. 

    %     %%% VERY CUMBERSOME APPROACH
    %     tmpout=[];
    %     for i=1:n_time_foldings
    %         tmpout=cat(1,tmpout,permute(Xs{3}(foldings_idx(i,:),:,:,:),[2 1 3 4])); %the ns_time_folding samples are concatenated along the second dimension
    %     end
    
    %%% less cumbersome way to fold iEEG
    disp('folding iEEG in time')
    tic
    tmp=Xs{3};
    tmp=permute(tmp(1:max(foldings_idx(:)),:,:,:),[1 5 2 3 4]);
    stmp=size(tmp);
    tmp=reshape(tmp,[ns_time_folding n_time_foldings stmp(3:end)]);
    tmp=permute(tmp,[2 1 4 5 3]);
    Xs{3}=tmp;
    disp('iEEG folded')
    toc
    
    %%% let's then fold the acoustic and semantic time-courses
    %%% and also add the lags (to be concatenated as additional models in
    %%% third dimension). 

    %     %%% VERY CUMBERSOME APPROACH
    %     disp('folding models in time, and adding lags')
    %     for acosem=1:2
    %         tmpout2=[];
    %         for ilag=1:n_lags
    %             tmpout1=[];
    %             for i=1:n_time_foldings
    %                 tmpout1=cat(1,tmpout1,permute(Xs{acosem}(foldings_idx(i,:)+ns_lags(ilag),:,:,:),[2 1 3 4])); %the ns_time_folding samples are concatenated along the second dimension
    %                 disp(num2str([i acosem ilag]))
    %             end
    %             tmpout2=cat(3,tmpout2,tmpout1);
    %             disp(num2str([acosem ilag]))
    %         end
    %         %         tmp1=tmpout2;
    %         Xs{acosem}=tmpout2;
    %     end
    
    %%% less cumbersome way to fold acoustics and semantics
    disp('folding acoustics and semantics in time')
    tic
    for acosem=1:2
        tmpout=[];
        for ilag=1:n_lags
            ns_thislag=ns_lags(ilag);
            tmp=Xs{acosem};
            tmp=permute(tmp([1:max(foldings_idx(:))]+ns_thislag,:,:,:),[1 5 2 3 4]); %#ok<NBRAK>
            stmp=size(tmp);
            tmp=reshape(tmp,[ns_time_folding n_time_foldings stmp(3:end)]);
            tmp=permute(tmp,[2 1 3 4 5]);
            stmp=size(tmp);
            tmp=reshape(tmp,[stmp(1) stmp(2)*stmp(3) stmp(4:5)]);
            if ilag==1
                stmp=size(tmp);
                stmp2=stmp;
                stmp2(3)=stmp2(3)*n_lags;
                tmpout=zeros(stmp2,'single'); %preallocate output matrix to speed up
            end
            tmpout(:,:,[1:size(tmp,3)]+(ilag-1)*size(tmp,3),:)=tmp; %#ok<NBRAK,SAGROW>
            %             disp(num2str([acosem ilag]))
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
Ds{1}=BLG_CosDistND(permute(Xs{1},[2 1 3 4]));
Ds{2}=BLG_CosDistND(permute(Xs{2},[2 1 3 4]));
Ds{3}=BLG_EucDistND(permute(Xs{3},[2 1 3 4]));
disp('Distances computed')
toc





















%% Acoustics|Semantics, single-subject plugin:
%% 1. shift in time the model, not the brain
%% 2. consider multiple model shifts as predictors (lags)
%% 3. do distances by considering nearby signal samples as dimensions
%% 4. do CV correlations here, not yet R^2_CV


rootmain='D:/!!Projects/DFI_iEEG/';
do_ft=1;
cd(rootmain)
run([rootmain,'code/project_init.m'])



Subjs={'subj-01'
    'subj-02'
    'subj-03'
    'subj-04'
    'subj-05'
    'subj-06'
    'subj-07'
    'subj-08'
    'subj-09'
    'subj-10'
    'subj-11'
    'subj-12'
    'subj-13'
    'subj-14'
    'subj-15'
    'subj-16'
    'subj-17'
    'subj-18'};
nsubjs=length(Subjs);

downsampfactor=20; %take one every downsampfactor samples of the 500 Hz timecourse
do_orth=1; %orthogonalize acoustics relative to semantics and vice versa
nswin=45; %chunk size in seconds
lags_ms=[0:20:300];%feature-brain-lags considered in ms (positive = feature before brain)

downsampfactor=50; %take one every downsampfactor samples of the 500 Hz timecourse
do_orth=1; %orthogonalize acoustics relative to semantics and vice versa
nswin=15; %chunk size in seconds
lags_ms=[0:20:100];%feature-brain-lags considered in ms (positive = feature before brain)
dodist=1; %do distance-based analyses


nlags=length(lags_ms);

audiofnout=[rootmain,'stimuli/AuditoryDimensions_Clean.mat'];
load(audiofnout,'X','fs','time','names')
t_aud=time;
dat_aud=X;

dat_aud_diff=dat_aud-[zeros(size(dat_aud(1,:)));dat_aud(1:end-1,:)]; %%% create temporal differential of audio features
dat_aud=cat(2,dat_aud,dat_aud_diff); %add differential to features

max_dur=9*60*1000;
idx_main=find(t_aud<max_dur);
max_ns=length(idx_main);
t_aud=t_aud(idx_main);
dat_aud=dat_aud(idx_main,:);
dat_aud(isinf(dat_aud))=0;


nam_aud=names';
nam_aud=cat(1,nam_aud,celfun(@(x)['d',x],nam_aud)); %add differential names
dat_aud=single(dat_aud);


t_aud=t_aud(1:2:end);
dat_aud=dat_aud(1:2:end,:);
disp(num2str(prctile(t_aud,[0 100]./1000./60)))
fs_aud=500;

dat_aud(isinf(dat_aud))=0;



lags_ns=unique(round(lags_ms.*fs_aud./1000));
idxout=find(t_aud>max(t_aud)-max(lags_ms));
t_aud(idxout)=[];
dat_aud(idxout,:)=[];
[aud_ns,aud_nfeat]=size(dat_aud);
win_ns=nswin*fs_aud;

dur_slice_ms=nswin*1000;
dur_slice_ns=dur_slice_ms/1000*fs_aud;

[slice_ini,slice_end]=ComputeSlices(aud_ns,dur_slice_ns);
check=slice_end(end-1:end)-slice_ini(end-1:end);
if std(check)~=0
    slice_ini(end)=[];
    slice_end(end)=[];
end

aud_idx=reshape(1:numel(dat_aud),size(dat_aud));
aud_idx_slice=[];

for i=1:size(dat_aud,2)
    tmp=aud_idx(:,i);
    tmp2=[];
    for j=1:length(slice_ini)
        tmp2=cat(2,tmp2,tmp(slice_ini(j):slice_end(j)));
    end
    aud_idx_slice=cat(3,aud_idx_slice,tmp2);
end

aud_idx_slice_resamp=[];
resampidx=[1:size(aud_idx_slice,1)]';
resampidx=resampidx(1:downsampfactor:end);
for i=0:downsampfactor-1
    tmp=aud_idx_slice(i+resampidx,:,:);
    aud_idx_slice_resamp=cat(4,aud_idx_slice_resamp,tmp);
end

aud_idx_slice_resamp=permute(aud_idx_slice_resamp,[4 1 2 3]);
size(aud_idx_slice_resamp)

%%% OK, let's estimate the acoustical distances!

d_aud=[];
disp('computing auditory dimension distances')
for i=1:length(lags_ns)
    tmp=dat_aud(aud_idx_slice_resamp+lags_ns(i));
    tmp=BLG_EucDistND(tmp);
    d_aud=cat(4,d_aud,tmp);
    clear tmp
end
disp('auditory dimension distances computed')
d_aud=permute(d_aud,[1 4 3 2]); %time pairs, lags, features, blocks, 
d_aud=real(d_aud);



if dodist==0
    adddistfn='nodist';
elseif dodist==1
    adddistfn='dist';
end
addfnam=['_slice_allsamp_',num2str(dur_slice_ms/1000),'s_win_',num2str(downsampfactor),'_downsamp_',num2str(nlags),'_nlags_',adddistfn];
if do_orth==0
    addfnam=['_corrs',addfnam];
elseif do_orth==1
    addfnam=['_pcorrs',addfnam];
end
fnsout=celfun(@(x) [rootmain,'participants/',x,'/stats/lagged',addfnam,'_aud_res_nlp.mat'],Subjs);


d_all={d_aud};
d_all=celfun(@(x) reshape(x,[size(x,1) size(x,2)*size(x,3) size(x,4)]),d_all); %collapse lags and features
clear d_aud d_nlp
nam_features=cat(1,nam_aud,nam_nlp);
nam_all={'acou' 'sem'};


d_all=celfun(@zscore,d_all);
%d_all contains acoustics distances.


fnsdat=celfun(@(x) [rootmain,'participants/',x,'/raw.fif'],Subjs);
out_dir=celfun(@(x) [rootmain,'participants/',x,'/stats/'],Subjs);
for i=1:length(out_dir)
    BLG_mkdir(out_dir{i})
end


for i=1:nsubjs
    warning off
    disp(['reading data for participant ',num2str([i nsubjs])])
    dat=ft_read_data(fnsdat{i})';
    
    head=ft_read_header(fnsdat{i});
    disp('data read')
    
    dat=single(dat);
    [nt_eeg,~]=size(dat);
    t_eeg=[0:nt_eeg-1]'./head.Fs.*1000;
    idxout=find(t_eeg<max(lags_ms));
    t_eeg(idxout)=[];
    dat(idxout,:)=[];
    [nt_eeg,nsens]=size(dat);
    
    disp('time: 0th and 100th percentile')
    disp(num2str(prctile(t_eeg,[0 100])./1000./60))
    disp(num2str(prctile(t_aud,[0 100])./1000./60))
    
    tmp_idx=reshape(1:numel(dat),size(dat));
    tmp_idx_slice=[];
    
    for ii=1:size(dat,2)
        tmp=tmp_idx(:,ii);
        tmp2=[];
        for j=1:length(slice_ini)
            tmp2=cat(2,tmp2,tmp(slice_ini(j):slice_end(j)));
        end
        tmp_idx_slice=cat(3,tmp_idx_slice,tmp2);
    end
    
    tmp_idx_slice_resamp=[];
    resampidx=[1:size(tmp_idx_slice,1)]'; %#ok<NBRAK>
    resampidx=resampidx(1:downsampfactor:end);
    for ii=0:downsampfactor-1
        tmp=tmp_idx_slice(ii+resampidx,:,:);
        tmp_idx_slice_resamp=cat(4,tmp_idx_slice_resamp,tmp);
    end
    
    tmp_idx_slice_resamp=permute(tmp_idx_slice_resamp,[4 1 2 3]);
    size(tmp_idx_slice_resamp)
    dat=dat(tmp_idx_slice_resamp);
    clear tmp_idx_slice_resamp tmp_idx_slice
    
    
    
    disp('computing iEEG distances')
    d_eeg=BLG_EucDistND(dat);
    d_eeg=single(d_eeg);
    d_eeg=zscore(d_eeg);
    disp('iEEG distances computed')
    d_eeg=permute(d_eeg,[1 4 2 3]);
    
    
    
    disp('computing GLMs and CV stats')
    
    disp('computing for acoustics')
    BetaAcou=BLG_GLM_ND(d_all{1},d_eeg,0,0); %training betas of acoustics-based model
    BetaAcouVal=BetaAcou;
    nb=size(d_eeg,3);
    for iblock=1:nb %average betas across non-test i.e., training blocks
        BetaAcouVal(:,:,iblock,:)=mean(BetaAcou(:,:,setxor(iblock,1:nb),:),3);
    end
    
    tmp=cat(2,d_all{1},ones(size(d_all{1}(:,1,:)))); %regression matrices of testing blocks
    pred=mtimesx(tmp,BetaAcouVal); %prediction on testing blocks based on betas in non-test (i.e., training) blocks
    cacou=(BLGmx_corr2(pred,d_eeg)); %cross-validated correlation between testing-block iEEG distance and training-based prediction
    
    disp('acoustics done')
    
    disp('computing for semantics')
    BetaNLP=BLG_GLM_ND(d_all{2},d_eeg,0,0); %training betas of semantics-based model
    BetaNLPVal=BetaNLP;
    nb=size(d_eeg,3);
    for iblock=1:nb
        BetaNLPVal(:,:,iblock,:)=mean(BetaNLPVal(:,:,setxor(iblock,1:nb),:),3);
    end
    
    tmp=cat(2,d_all{2},ones(size(d_all{2}(:,1,:))));%regression matrices of testing blocks
    pred=mtimesx(tmp,BetaNLPVal);%prediction on testing blocks based on betas in non-test (i.e., training) blocks
    cnlp=(BLGmx_corr2(pred,d_eeg));%cross-validated correlation between testing-block iEEG distance and training-based prediction
    
    disp('semantics done')
    
    %output CV correlations
    % matrix of CV correlations c of size [[aco,sem] nwindows nsensors]
    c=squeeze(cat(1,cacou,cnlp));
    
    clear pred d_eeg
    disp(['all done: ',num2str([i nsubjs])])
    save(fnsout{i},'c','lags_ns','fs_aud','lags_ms','nam_all','nam_features')
    
end

