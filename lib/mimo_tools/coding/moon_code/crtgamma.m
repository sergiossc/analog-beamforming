function gamma = crtgamma(m)
% function gamma = crtgamma(m)
%
% Compute the gammas for the CRT

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

r = length(m);
mp = prod(m);
for i=1:r
  [g,b,y1] = gcdint2(mp/m(i),m(i));
  gamma(i) = mp/m(i)*b;
end