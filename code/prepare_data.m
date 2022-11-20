
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

