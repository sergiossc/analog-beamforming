function v = voln(n)
% function v = voln(n)
% Compute the volume of an n-dimensional sphere
%

% Todd Moon, March 20, 2004
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

if(~mod(n,2)) % even
  v = pi^(n/2)/factorial(n/2);
else
  v = 2^n*pi^((n-1)/2)*factorial((n-1)/2)/factorial(n);
end