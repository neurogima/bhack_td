function Y=BLG_diagND(X)
% function Y=BLG_diagND(X)
% given in input a distance matrix X of size [Nobjs Nobjs N3 N4... Nn]
% where NobjsxNobjs is symmetric
% returns in output a matrix Y of size [Nobjs N3 N4,... Nn]
% where each column contains the diagonal elements of each square page
% Bruno L. Giordano
% brungio@gmail.com
% November 2017


%% debug and sanity checks
% X=rand([50 50 10 2 4]);

sX=size(X);
if sX(1)~=sX(2)
    error('X pages must be square matrices')
end

nobjs=sX(1);
% npairs=nobjs*(nobjs-1)/2;
sXrest=sX(3:end);
if isempty(sXrest)
    sXrest=1;
end
nXrest=prod(sXrest);
X=reshape(X,[nobjs nobjs nXrest]);

%%% the loop is actually the quickest!
% %%% debug and sanity checks
% tic
Y=zeros(nobjs,nXrest);
for i=1:nXrest
    Y(:,i)=diag(X(:,:,i));
end
% Ycheck=Y;
%%% debug and sanity checks
% toc


% %%% debug and sanity checks
% % tic
% % Y=zeros(nobjs,nXrest);
% idx=reshape(1:nobjs^2,[nobjs nobjs]);
% idx=diag(idx);
% idx=repmat(idx,[1 nXrest]);
% addidx=[0:prod(sXrest)-1]*nobjs*nobjs;
% idx=bsxfun(@plus,idx,addidx); %indices to what we want in the right order
% Y=X(idx);
% %%% debug and sanity checks
% % toc
% isequal(Ycheck,Y)

Y=reshape(Y,[nobjs sXrest]);
% 
% %%
% nperms=1000;
% 
% 
% x=rand([100 2 nperms]);
% D=BLG_AbsDist(x);
% r=BLGmx_corr(x);
% rD=BLGmx_corr(D);
% 
% r=BLG_vectorformND(r)';
% rD=BLG_vectorformND(rD)';
% 
% close all
% hist([r rD],100)
% 
% 
% 















