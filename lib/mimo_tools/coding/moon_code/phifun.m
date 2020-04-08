function phi = phifun(x,sigma)
% function phi = phifun(x,sigma)
%
% Todd K. Moon, Sept. 30, 2003

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

s22 = 1/(2*sigma*sigma);
c = 1/sqrt(8*pi*sigma*sigma);

phi = c*(exp(-(x-1).^2*s22) + exp(-(x+1).^2*s22));
