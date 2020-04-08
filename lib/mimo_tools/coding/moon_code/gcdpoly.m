function [g,s,t] = gcdpoly(b,c,thresh)
% function [g,s,t] = gcdpoly(b,c)
% Compute the GCD g = (b,c) using the Euclidean algorithm
% and return s,t such that bs+ct = g, where b and c are polynomials
% with real coefficients
%
% thresh = (optional) threhold argument used to truncate small remainders

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

rm2 = b; rm1 = c;
sm2 = 1; sm1 = 0;
tm2 = 0; tm1 = 1;
while(any(rm1))
  [q,tr] = polydiv(rm2,rm1);
  if(nargin==3)	tr(find(abs(tr) < thresh)) = 0; end; % truncate small
  ts = polysub(sm2,polymult(q,sm1));
  tt = polysub(tm2,polymult(q,tm1));
  rm2 = rm1; sm2 = sm1;  tm2 = tm1; 
  rm1 = tr;  sm1 = ts;  tm1 = tt;
end
lc = rm2(1);					% make monic
g = rm2/lc;
s = sm2/lc;
t = tm2/lc;
