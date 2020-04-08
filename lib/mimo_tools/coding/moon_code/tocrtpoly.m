function y = tocrt(f,m)
% function y = tocrt(f,m)
%
% Compute the representation of the polynomial f using the
% using the Chinese Remainder Theorem (CRT) with
% moduli m = {m1,m2,...,mr}.  It is assumed (without checking)
% that the moduli are relatively prime.
% m is passed in as a cell array containing polynomial vectors
% and y is returned as a cell array containing polynomial vectors

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[n,r] = size(m);
for i=1:r
  [q,y{i}] = polydiv(f,m{i});
end

