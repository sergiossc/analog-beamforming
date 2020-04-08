function [lu,indx] = lu2(A)
% function [lu,indx] = lu2(A)
%
% Compute the lu factorization of a binary matrix A
%
% lu = matrix contining L and U factors
% indx = index of pivot permutations
%
% If A is not square, it is assumed that it is tall, 
% and the LU is computed in such a way that PA has an
% inverse in the top colsxcols matrix

% Todd K. Moon
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[m,n] = size(A);
if(m < n)
  error('Matrix should be square or tall');
end
indx = 1:m;
for k=1:n-1
  % for pivoting determine the largest element in this column
  % midx = index of largest
  [p,midx] = max(abs(A(k:m,k)));
  % the previous index was out of k:n; adjust it so it is
  % index on 1:n
  midx = midx+k-1;
  % interchange the mth and (midx)th rows
  dum = A(k,1:n);  A(k,1:n) = A(midx,1:n);  A(midx,1:n) = dum;
  % record on which row the kth row was swapped
  dum1 = indx(k);
  indx(k) = indx(midx);
  indx(midx) = dum1;
  if(A(k,k) == 0)
    error('Linearly dependent columns');
  else
    for j=k+1:m
      mult = A(j,k);
      % do the row operation
      A(j,k:n) = mod(A(j,k:n) - mult*A(k,k:n),2);
      % store the multiplier elment in the lower triangle
      A(j,k) = mult;
    end
  end
end
lu = A;