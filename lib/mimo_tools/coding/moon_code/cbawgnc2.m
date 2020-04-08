function c = cbawgnc2(Ecsigma2)
% function c = cbawgnc2(Ecsigma2)
% Compute the capacity for a BAWGNC at a coded signal to noise ratio
% of Ec/sigma^2

% Assume that Ec = 1, so sigma^2 = 1/Ecsigma2

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

sigma2 = 1/Ecsigma2;
sigma = sqrt(sigma2);

LOLIM = -30*sigma;
UPLIM = 30*sigma;
e = exp(1);

c = -quad(@philog,LOLIM,UPLIM,[],[],sigma) - 0.5*log2(2*pi*e*sigma2);
