% function out = celfun(x,y)
% same as cellfun, with 'un' = 0
% Bruno L. Giordano
% brungio@gmail.com
% November 2018
function out = celfun(x,y)

out=cellfun(x,y,'un',0);