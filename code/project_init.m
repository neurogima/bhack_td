clc
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Preparing the environment')

cd(rootmain)
disp(['rootmain = ',rootmain])

restoredefaultpath
disp('Default path restored')
% 
% if ~exist('do_ft','var') || do_ft
%     ftpath='D:\Matlab_stuff\fieldtrip-20200828\';
%     addpath(ftpath)
%     ft_defaults
%     disp('FT added to path')
% end
% 
% if ~exist('do_spm','var') || do_spm
%     spmpath='D:\Matlab_stuff\spm12\';
%     addpath(spmpath)
%     disp('SPM12 added to path')
%     
%     addpath(genpath('D:\!!Projects\ToolboxRepo\xjview97'));
%     disp('xjview added to path')
% end
% 

addpath(genpath([rootmain,'code/']))
disp('Project code repo added to path')

disp('WARNING: random seeds not initialized')
disp('(required for e.g. proper permutation)')

disp('All ready, wise master')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

