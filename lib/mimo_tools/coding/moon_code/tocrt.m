function y = tocrt(x,m)
% function y = tocrt(x,m)
%
% Compute the representation of the scalar x using the
% using the Chinese Remainder Theorem (CRT) with
% moduli m = [m1,m2,...,mr].  It is assumed (without checking)
% that the moduli are relatively prime

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

y = mod(x,m);

