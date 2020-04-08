function [x,gamma] = fromcrt(y,m,gammain)
% function x = fromcrt(y,m)
%
% Given a sequence [y1,y2,...,y2] that is a representation 
% of an integer x in CRT form, convert back to x.
% m = [m1,m2,...,mr]
% gammain (optional) = set of gamma factors.  If not passed in, 
% gamma is computed
%

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

r = length(y);
mp = prod(m);
if(nargin==2)
  x = 0;
  for i=1:r
    [g,b,y1] = gcdint2(mp/m(i),m(i));
    gamma(i) = mp/m(i)*b;
    x = x + gamma(i)*y(i);
  end
else
  x = gammain*y';
  gamma = gammain;
end
x = mod(x,mp);