function philo = philog(x,sigma)
% function phi = philog(x,sigma)
%
% Todd K. Moon, Sept. 30, 2003

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

phi = phifun(x,sigma);
philo = phi.*log2(phi);
