function r=BLGmx_corr2(x,y)
% why is Matlab scared of pages?
%
% function r=BLGmx_corr2(x,y)
% replaces Matlab's corr(x,y)
% requires mtimesx
%
% tip: reshape input matrices so that prod(3:,...ndims) matrix dimensions
% are all in third dimension of input. Then apply appropriate inverse
% reshape to r, if required.
%
% Bruno L. Giordano brungio@gmail.com
% INP, Glasgow
% INT, CNRS, Marseille
% 2017: first version.
% 2018: accepts x and y; quicker bsxfun approach from Matlab.

if nargin<2
    corrXX=1;
end

sx=size(x);
n=sx(1);
nx=sx(2);
nxdims=length(sx);

x = bsxfun(@minus,x,sum(x,1)/n);  % demean
if nargin < 2
    r = mtimesx(x,'t',x,'n');
    d = sqrt(BLG_diagND_noSQ(r));
    r = bsxfun(@rdivide,r,d);
    r = bsxfun(@rdivide,r,permute(d,[2 1 3:nxdims]));
else
    ny=size(y,2);
    y = bsxfun(@minus,y,sum(y,1)/n);  % demean
    r = mtimesx(x,'t',y,'n'); 
    dx = sqrt(sum(abs(x).^2, 1));
    dy = sqrt(sum(abs(y).^2, 1));
    r = bsxfun(@rdivide,r,permute(dx,[2 1 3:nxdims]));
    r = bsxfun(@rdivide,r,dy);
end
% BLG: potentially do this, takes time and diagonal signs are more often
% than not uninteresting.
% % Limit off-diag elements to [-1,1], and put exact ones on the diagonal for autocorrelation.
% t = find(abs(r) > 1); r(t) = r(t)./abs(r(t)); % preserves NaNs
% if corrXX
%     r(1:nv+1:end,:,:) = sign(BLG_diagND_noSQ(r)); % preserves NaNs
% end

function Y=BLG_diagND_noSQ(X)
% why is Matlab scared of pages?
%
% function Y=BLG_diagND_noSQ(X)
% given in input a distance matrix X of size [Nobjs Nobjs N3 N4... Nn]
% where NobjsxNobjs is symmetric
% returns in output a matrix Y of size [Nobjs 1 N3 N4,... Nn]
% where each column contains the diagonal elements of each square page
%
% Bruno L. Giordano brungio@gmail.com
% INP, Glasgow
% INT, CNRS, Marseille
% _no squeeze(Y) edit:September 2018

sX=size(X);
if sX(1)~=sX(2)
    error('X pages must be square matrices')
end

nobjs=sX(1);
sXrest=sX(3:end);
if isempty(sXrest)
    sXrest=1;
end
nXrest=prod(sXrest);
X=reshape(X,[nobjs nobjs nXrest]);

%%% the loop is actually the quickest!
Y=zeros(nobjs,nXrest);
for i=1:nXrest
    Y(:,i)=diag(X(:,:,i));
end
Y=reshape(Y,[nobjs 1 sXrest]);