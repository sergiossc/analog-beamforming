function x = numtonum(xsi,SkINTSi,Q)
% given an array, extract the positions indicated in the array SkINTSi
% and convert to a base Q number

x = 0;
for i=length(SkINTSi):-1:1
   x = Q*x + xsi(SkINTSi(i));
end
