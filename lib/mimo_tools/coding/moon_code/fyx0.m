function f = fyx0(x,ysigma2)
% compute f(y|x)
% the value of x is a 0 (meaning a -1 is sent)
% or a 1 (meaning a 1 is sent)
% The argument ysigma2 is an array of "auxiliary" values, 
% ysigma2(1) = y value,
% ysigma2(2) = sigma^2
% ysigma2(3) = a

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

y = ysigma2(1);
sigma2 = ysigma2(2);
a = ysigma2(3);

if(x==0) 
  f = exp(-1/(2*sigma2)*(y+a)^2);
elseif(x==1) 
  f = exp(-1/(2*sigma2)*(y-a)^2);
end