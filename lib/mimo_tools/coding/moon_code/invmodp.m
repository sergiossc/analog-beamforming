function binv = invmodp(b, p)
% function binv = invmodp(b, p)
% 
% Compute the inverse of b modulo p

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

if(b == 0)
  error('Error: Illegal division by zero');
end
[g,x,y] = gcd(p,b);
if(g ~= 1)
  error('Error: attempting to divide by noninvertible number');
end
binv =  mod(y,p);

