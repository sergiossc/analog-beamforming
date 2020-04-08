function Ainv = inv2(A)
%function Ainv = inv2(A)
% Compute the inverse of a binary matrix A

% Todd Moon
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[m,n] = size(A);
if(m < n)
  error('Matrix should be square or tall');
end
[lu,ind] = lu2(A);
ind
b = zeros(m,1);
Ainv = [];
for i=1:n
  if(i>1) b(i-1) = 0; end;
  b(i) = 1;
  Ainv = [Ainv lubksub2(lu,ind,b)];
end