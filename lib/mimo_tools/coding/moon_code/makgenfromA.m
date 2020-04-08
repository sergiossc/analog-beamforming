function [G,success,A] = makegenfromA(A)
% function [G,success,A] = makegenfromA(A)
% given an MxN parity check matrix A (not necessarily systematic)
% with N>M, compute a systematic generator matrix G

% Todd K. Moon
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[M,N] = size(A);
K = N-M;
[success,Ainv,Ac,pidx] = gaussj2(A);
A = A(:,pidx);

% mod(Ainv*A,2)
G = [Ac; eye(K)];
