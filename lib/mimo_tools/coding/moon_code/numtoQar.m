function x = numtoQar(xint,Q,Neval,varlist,Nvar)
% given a number xint which represents a base-Q representation
% of a number with Neval digits, find the base-Q representation.  The 
% varlist indicates which positions the
% digits should go into, out of a possible set of Nvar.
% The least significant bits are placed in the variable positions
% listed FIRST in varlist

x = zeros(1,Nvar);
for i=1:Neval
  d = mod(xint,Q);
  xint = floor(xint/Q);
  x(varlist(i)) = d;
end
