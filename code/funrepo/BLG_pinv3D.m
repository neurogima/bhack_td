function Y=BLG_pinv3D(X,dopar)
%%% Y=BLG_pinv3D(X)
%%% faster pseudo-inverse of each page of X
%%% through mtimesx and elimination of
%%% unnecessary loops
%%% speed boost only for square pages now
%%%
%%% Bruno L. Giordano
%%% brungio@gmail.com
%%% March 2018

if nargin<1
    dopar=0; %no parallel computations for svd
end

[x1,x2,NP]=size(X);
if x1~=x2
    % tic
    Y=permute(X,[2 1 3]);
    for i=1:NP
        Y(:,:,i)=pinv(X(:,:,i));
    end
    % toc
else
    % tic
    U=X;
    S=X;
    V=X;
    
    %p=gcp('nocreate');
    if dopar%~isempty(p)
        parfor i=1:NP
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i),'econ');
        end
    else
        for i=1:NP
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i),'econ');
        end
    end
    
    
    
    s=permute(BLG_diagND(S),[3 1 2]);
    tol=max(size(X(:,:,1))).*eps(max(abs(s),[],2));
    r1 = sum(bsxfun(@minus,s,tol)>0)+1;
    ur1=unique(r1);
    nur1=length(ur1);
    
    s = 1./s;
    if nur1==1
        V(:,ur1:end,:) = [];
        U(:,ur1:end,:) = [];
        s(:,ur1:end,:) = [];
        Y=mtimesx(bsxfun(@times,V,s),'n',U,'t');
    else
        Y=X;
        for i=1:nur1
            thisur1=ur1(i);
            idx=r1==thisur1;
            thisV=V(:,1:thisur1-1,idx);
            thisU=U(:,1:thisur1-1,idx);
            thiss=s(:,1:thisur1-1,idx);
            Y(:,:,idx)=mtimesx(bsxfun(@times,thisV,thiss),'n',thisU,'t');
        end
    end
    %toc
end

%%

%
% % % normS=arrayfun(@(idx) norm(S(:,:,idx),inf),1:NP,'un',false);
% % tic
% % [U,S,V]=arrayfun(@(idx) svd(X(:,:,idx)),1:NP,'un',false);
% % toc
% %
% % tic
% % Y=arrayfun(@(idx) pinv(X(:,:,idx)),1:NP,'un',false);
% % toc
% %%
%
% function X = pinv(A,tol)
% %PINV   Pseudoinverse.
% %   X = PINV(A) produces a matrix X of the same dimensions
% %   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
% %   are Hermitian. The computation is based on SVD(A) and any
% %   singular values less than a tolerance are treated as zero.
% %   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS(class(A)).
% %
% %   PINV(A,TOL) uses the tolerance TOL instead of the default.
% %
% %   Class support for input A:
% %      float: double, single
% %
% %   See also RANK.
%
% %   Copyright 1984-2013 The MathWorks, Inc.
%
% [U,S,V] = svd(A,'econ');
% s = diag(S);
% if nargin < 2
%     tol = max(size(A)) * eps(norm(s,inf));
% end
% r1 = sum(s > tol)+1;
% V(:,r1:end) = [];
% U(:,r1:end) = [];
% s(r1:end) = [];
% s = 1./s(:);
% X = bsxfun(@times,V,s.')*U';
