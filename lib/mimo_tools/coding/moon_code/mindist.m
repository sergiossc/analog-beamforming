function dmin = mindist(G)
% function dmin = mindist(G)
% G is an NxK matrix
%
% Compute the minimum distance of a binary code
% by exhaustive checking.  (Therefore, only good for codes up to 
% moderate sizes

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[N,K] = size(G);
dmin = N;
m = zeros(K,1);
for i=1:(2^K-1)
  m = c2b(i,K)';
  c = mod(G*m,2);
  d = sum(c);
  if(d < dmin)
	dmin = d;
	mmin = m;
	cmin = c;
  end
end

dmin
mmin
cmin

  
  