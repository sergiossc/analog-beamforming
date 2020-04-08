function c = d2b(n,m)
% function c = d2b(n,m)
% convert n to an m-bit binary representation

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

mask = 1;
c = zeros(1,m);
for i=1:m
  if(bitand(mask,n))
	c(i) = 1;
  end
  mask = mask*2;
end
	
