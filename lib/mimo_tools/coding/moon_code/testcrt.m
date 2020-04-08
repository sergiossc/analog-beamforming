% An example of CRT calculations

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

m = [4 27 25];
gamma = crtgamma(m);  % precompute the gamma values

a = [0 2 3];
x = fromcrt(a,m);

