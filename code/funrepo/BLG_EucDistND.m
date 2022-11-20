function D=BLG_EucDistND(X,keepsquare)
% function D=BLG_EucDistND(X)
% - X = data matrix of dimensions:
%   1. NX1 = n elements for distance computation (e.g., voxels in fMRI) - EEG
%      channels (or time points)
%   2. Nconditions (or stimuli or etc): the distances are computed between
%       what in this dimension
%   3. NDists = number of distances
%
%   D = distance matrix of dimensions:
%   1. NPairs = ncondition pairs (ordered as expected by squareform)
%   2. NDists
% - keepsquare = 1 if the output is instead square distance matrices
%
%   Accepts input matrixes X of ndimensions > 3. In this case, D will have
%   more than two dimensions, e.g., 
%   if size(X) = [NX1 NConditions NY3 NY4]
%   then size(D) = [NConditionPairs NY3 NY4] if keepsq=0
%                  [NConditions NConditions NY3 NY4] if keepsq=1
%
%   Bruno L. Giordano
%   December 2017
%   brungio@gmail.com
%%
% clc
% X=rand([100 15 100 5]);
if nargin<2
    keepsquare=0;
end

sX=size(X);
NX1=sX(1);
NConds=sX(2);
if length(sX)==2
    sRest=1;
else
    sRest=sX(3:end);
    X=reshape(X,[NX1 NConds prod(sRest)]);
end

% tic

C=mtimesx(X,'t',X,'n');
Cdiag=BLG_diagND(C);
Cdiag=permute(Cdiag,[1 3 2]);
D=bsxfun(@plus,Cdiag,permute(Cdiag,[2 1 3]));
D=D-2*C;
D=sqrt(D);

if ~keepsquare
    D=BLG_vectorformND(D);
    D=reshape(D,[size(D,1) sRest]);
else
    D=reshape(D,[size(D,1) size(D,2) sRest]);
end



% toc
% 
% tic
% d=zeros(size(D));
% for i=1:prod(sRest)
%     d(:,i)=pdist(X(:,:,i)')';
% end
% toc
% 
% close all
% plot(D(:),d(:),'.')



