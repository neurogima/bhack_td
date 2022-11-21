function [Beta,Pred,RSQ,Xinv]=BLG_GLM_ND(X,Y,dozscore,parloopflag)
% faster mtimesx-based GLMs with support for pages.
%
% X = predictors matrix of size [ndata npredictors dim3x, dim4x,... dimnx]
% Y = dependent variable matrix of size [ndata npredictors dim3y, dim4y,... dimny]
% dozscore = zscore X and Y
% parloopflag = compute inverse for different pages within parfor loop
%
% outputs:
% Beta = GLM coefficients
% Pred = prediction
% RSQ = R squared
% Xinv = inverse
%
%
% Handy features:
% - input X doesn't need to include an intercept term.
% - singleton dimensions>2 allowed and automatically expanded to match
%   non-singleton dimensions in the other variable
%  e.g.:
% X=rand([100,10,1,4]);
% Y=rand([100,1,10]);
% [Beta,Pred,RSQ,Xinv]=BLG_GLM_ND(X,Y)
% produces:
% Betas of size 11 1 10 4
%
% Bruno L. Giordano
% brungio@gmail.com
% November 2018
% June 2021: - implement automatic expansion of singleton dimensions>2 of X/Y
%              for corresponding non-singleton dimensions>2 in Y/X
%            - input matrices can have more than 3 dimensions


if nargin<3
    dozscore=0; %default
end

if nargin<4
    parloopflag=0; %default
end

%debug/devel
% dozscore=0;
% parloopflag=0;
% X=rand([100,10,1,4]);
% Y=rand([100,1,10]);

sX=size(X);
sY=size(Y);

if sX(1)~=sY(1)
    error('n datapoints for X and Y differs')
end

if sY(2)~=1
    error('only one predicted Y variable allowed')
end

nXdims=length(sX);
nYdims=length(sY);
if nXdims<nYdims
    sX=[sX ones(1,nYdims-nXdims)];
elseif nYdims<nXdims
    sY=[sY ones(1,nXdims-nYdims)];
end

if length(sX)==2
    sX=[sX 1];
end
if length(sY)==2
    sY=[sY 1];
end

% if any((sX(3:end)-sY(3:end))~=0)
%     warning('expanding variables along singleton pages')
% end

nXdims=length(sX);
nYdims=length(sY);
Xrep=ones(size(sX));
Yrep=ones(size(sY));

for i=3:nXdims
    if sX(i)==1 && sY(i)~=1
        Xrep(i)=sY(i);
    elseif sX(i)~=1 && sY(i)==1
        Yrep(i)=sX(i);
    elseif sX(i)~=1 && sY(i)~=1 && sX(i)~=sY(i)
        error('check consistency of X and Y pages')
    end
end


if any(Xrep~=1)
    redX=X;
end
X=repmat(X,Xrep);
Y=repmat(Y,Yrep);
sX=size(X);
sY=size(Y);

X=reshape(X,[sX(1:2),prod(sX(3:end))]);
Y=reshape(Y,[sY(1:2),prod(sY(3:end))]);

%%% sanity checks
% X=rand(100,10);
% Y=rand(100,1);
% [Beta,Pred,RSQ]=BLG_GLM_ND(X,Y);
% Xint=cat(2,ones(size(X(:,1,:))),bsxfun(@minus,X,mean(X)));
%
% [B,BINT,R,RINT,STATS] = regress(Y,Xint);
% corr(B,Beta)
% corr(Xint*B,Pred)
% disp(num2str([RSQ STATS(1)]))


fdemean=@(x) bsxfun(@minus,x,mean(x));
fzscore=@(x) zscore(x);
% fpinv=@(x) BLG_pinv3D(x,1);
feye=@(x,npages) repmat(eye(size(x,1)),[1 1 npages]);
% fint=@(s) ones(s); %add intercept
fresmatr=@(x,xpinv) feye(x,size(x,3))-mtimesx(x,xpinv); %residual forming matrix
fbeta=@(xinv,y) mtimesx(xinv,'n',y,'n'); %betas
fregrzsc=@(x) cat(2,ones(size(x(:,1,:))),fzscore(x)); %zscore and add intercept
fregrdem=@(x) cat(2,ones(size(x(:,1,:))),fdemean(x)); %demean and add intercept
fpred=@(x,beta) mtimesx(x,beta); %predictions
fres=@(xresmatr,y) mtimesx(xresmatr,y); %residuals
% frsq=@(y,ypred) BLG_vectorformND(BLGmx_corr(cat(2,y,ypred))).^2;
frsq=@(y,ypred) BLGmx_corr2(y,ypred).^2;

if dozscore
    X=fregrzsc(X);
    Y=fzscore(Y);
else % else just demean predictors
    X=fregrdem(X);
end
if any(Xrep~=1)
    if dozscore
        redX=fregrzsc(redX);
    else % else just demean predictors
        redX=fregrdem(redX);
    end
    sredx=size(redX);
    if length(sredx)==2
        sredx=[sredx 1];
    end
    redX=reshape(redX,[sredx(1:2) prod(sredx(3:end))]);
    % tmpXinv=zeros([sredx([2 1]) prod(sredx(3:end))]);
    tmpXinv=BLG_pinv3D(redX,parloopflag);
    tmpXinv=reshape(tmpXinv,[sredx([2 1]) sredx(3:end)]);
    tmpXinv=repmat(tmpXinv,Xrep);
    stmp=size(tmpXinv);
    Xinv=reshape(tmpXinv,[stmp(1:2) prod(stmp(3:end))]);
else
    Xinv=BLG_pinv3D(X,parloopflag);
end
Beta=fbeta(Xinv,Y);

if nargout>1
    Pred=fpred(X,Beta);
end
if nargout>2
    RSQ=frsq(Y,Pred);
end

if length(sX)>2
    Beta=reshape(Beta,[size(Beta,1) size(Beta,2) sX(3:end)]);
    if nargout>1
        Pred=reshape(Pred,[size(Pred,1) size(Pred,2) sX(3:end)]);
    end
    if nargout>2
        RSQ=reshape(RSQ,[size(RSQ,1) size(RSQ,2) sX(3:end)]);
    end
end