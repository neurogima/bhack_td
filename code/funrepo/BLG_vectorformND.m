function Y=BLG_vectorformND(X)
% function Y=BLG_vectorformND(X)
% given in input a distance matrix X of size [Nobjs Nobjs N3 N4... Nn]
% where NobjsxNobjs is symmetric
% returns in output a matrix Y of size [Nobjectpairs N3 N4,... Nn]
% where each column sorts between-object distances as expected by
% squareform, vectorform etc.
% Bruno L. Giordano
% brungio@gmail.com
% November 2017


%%% debug and sanity checks
% X=rand([36 36 1000]);
% 
% clc
% tic
% Ycheck=[];
% for i=1:size(X,3)
%     tmp=vectorform(X(:,:,i))';
%     X(:,:,i)=squareform(tmp);
%     Ycheck=cat(2,Ycheck,tmp);
% end
% toc

sX=size(X);
if sX(1)~=sX(2)
    error('X pages must be square matrices')
end
nobjs=sX(1);
npairs=nobjs*(nobjs-1)/2;
sXrest=sX(3:end);
if isempty(sXrest)
    sXrest=1;
end
nXrest=prod(sXrest);
X=reshape(X,[nobjs nobjs nXrest]);

%%% debug and sanity checks
% tic
% Y=zeros(npairs,nXrest);
% idx=reshape(1:nobjs^2,[nobjs nobjs]);
% idx=idx(tril(true(nobjs),-1));
% for i=1:nXrest
%     tmp=X(:,:,i);
%     Y(:,i)=tmp(idx);
%     %     X(:,:,i)=squareform(tmp);
% end
% toc

%%% debug and sanity checks
% tic
Y=zeros(npairs,nXrest);
idx=reshape(1:nobjs^2,[nobjs nobjs]);
idx=idx(tril(true(nobjs),-1));
idx=repmat(idx,[1 nXrest]);
addidx=[0:prod(sXrest)-1]*nobjs*nobjs;
idx=bsxfun(@plus,idx,addidx); %indices to what we want in the right order
Y=X(idx);
Y=reshape(Y,[npairs sXrest]);
%%% debug and sanity checks
% toc
% isequal(Ycheck,Y)



end



function Z = vectorform(Y)
[~, n] = size(Y);
Z = Y(tril(true(n),-1));
Z = Z(:)';% force to a row vector, even if empty
end

