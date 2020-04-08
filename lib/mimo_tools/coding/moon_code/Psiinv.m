function p = Psiinv(xlist,mlist2,psi2)
% function p = Psiinv(x,mlist2,psi2)
%
% Compute the Psi^{-1} function used in density evolution
% x = value
% mlist2, psi2 -- list of points to interpolate (see densev1)

% For large values, use the approximation from
% [Chung,Richardson,Urbanke2001]

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

p = [];

for x = xlist
  if(x==0)
	p1 = 0;
  end
  s = sign(x);
  x = abs(x);
  if(x > psi2(end)) 					% use function approximation
    p1 = ((log(1-x)-.0218)/(-.4527))^(1/.86);
  else
   p1 = interp1(psi2,mlist2,x);   % interpolate from the given data
  end
  p1 = p1*s; 								% correct the sign
  p = [p p1];
end
