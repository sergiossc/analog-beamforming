function gamma = crtgammapoly(m)
% function gamma = crtgammapoly(m)
%
% Compute the gammas for the CRT for polynomials

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

r = length(m);
mp = 1;
for i=1:r
  mp = polymult(mp,m{i});
end
for i=1:r
  [q,rm] = polydiv(mp,m{i});
  [g,b,y1] = gcdpoly(q,m{i});
  gamma{i} = polymult(q,b);
end