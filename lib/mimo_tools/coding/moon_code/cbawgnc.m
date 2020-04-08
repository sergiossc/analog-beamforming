function C = cbawgnc(sigmain)
% function C = cbawgnc(sigmain)

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

C = zeros(1,length(sigmain));
e = exp(1);
i = 0;
for sigma = sigmain
  i = i+1;
  LOLIM = -30*sigma;
  UPLIM = 30*sigma;
  C(i) = -quad(@philog,LOLIM,UPLIM,[],[],sigma) - 0.5*log2(2*pi*e*sigma*sigma);
end