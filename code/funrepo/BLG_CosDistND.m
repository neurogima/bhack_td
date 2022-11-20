function D=BLG_CosDistND(X)
% function D=BLG_CosDistND(X)
% - X = data matrix of dimensions:
%   1. NX1 = n elements for distance computation (e.g., voxels in fMRI) - EEG
%      channels (or time points)
%   2. Nconditions (or stimuli or etc): the distances are computed between
%       what in this dimension
%   3. NDists = number of distances
%
%   D = cosine distance matrix of dimensions:
%   1. NPairs = ncondition pairs (ordered as expected by squareform)
%   2. NDists
%
%   Accepts input matrixes X of ndimensions > 3. In this case, D will have
%   more than two dimensions, e.g., 
%   if size(X) = [NX1 NPairs NY3 NY4]
%   then size(D) = [NPairs NY3 NY4]
%
%   Bruno L. Giordano
%   July 2020
%   brungio@gmail.com
%%
% clc
% X=rand([100 15 100 5]);
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
X=bsxfun(@minus,X,mean(X));
SS1=sum(X.^2);
X=bsxfun(@times,X,1./sqrt(SS1));
D=1-mtimesx(X,'t',X,'n');
D=BLG_vectorformND(D);
% toc

% tic
% d=zeros(size(D));
% for i=1:prod(sRest)
%     d(:,i)=pdist(X(:,:,i)','cosine')';
% end
% toc
% 
% close all
% plot(D(:),d(:),'.')

D=reshape(D,[size(D,1) sRest]);


