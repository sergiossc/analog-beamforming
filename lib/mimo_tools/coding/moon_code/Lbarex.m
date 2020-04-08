%
% test the average performance of a GS(m) code
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

n = 31;
q = 32;
k = 15;
Mlist = [0,3,21];

for m=Mlist
  tm = computetm(n,k,m);
  Lm = computeLm(n,k,m);
  L = computeLbar(n,k,q,tm);
  fprintf(1,'m=%d  t=%d  Lm=%d  L=%g\n',m,tm,Lm,L);
end
  