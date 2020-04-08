function x = lubksub2(lu,indx,b)
% Solve Ax = b, where A has been factored into (lu,indx) using lu2

% Todd K. Moon
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[m,n] = size(lu);
if(m < n)
  error('Matrix should be square or tall');
end
b = b(:);
b =  b(indx);
% Do the backsubstitution to solve Ly = Pb, where P is determined by indx
x = zeros(n,1);
for i=1:n
  sum = b(i);
  j = 1:i-1;
  sum = sum - lu(i,j)*b(j);
  b(i) = mod(sum,2);
end
% Solve Ux = y
for i=n:-1:1
  sum = b(i);
  j = i+1:n;
  sum = sum - lu(i,j)*x(j);
  x(i) = mod(sum,2);
end