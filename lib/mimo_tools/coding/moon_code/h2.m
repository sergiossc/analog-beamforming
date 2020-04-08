function h = h2(p)
% function h = h2(p)
% Compute the binary entropy function

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

i1 = find(p==0 | p==1);
i2 = find(p ~=0 & p ~= 1);
h(i1) = zeros(size(i1));
h(i2) = -(p(i2) .* log2(p(i2)) + (1-p(i2)) .* log2(1-p(i2)));
