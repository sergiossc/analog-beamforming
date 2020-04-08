function  [SNRdB,sigma] = caprate(channel,r)
% function [SNRdB,sigma] = caprate(channel,r);
% 
% determine the capacity of either an awgn channel or
% a binary awgn channel.
% 
% This is done by solving for the sigma for which capacity(sigma) = rate
% where capacity(.) is the appropriate capacity function
%
% channel is either 'cbawgnc' or 'cawgnc'
% e.g.  
% caprate('cbawgnc',1/2)

% Todd K. Moon, May 7, 2004
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

% See, e.g., Richardson, Urbanke, p. 118, or p. 146)
  y = r;
  fun = inline(sprintf('%s(sigma)-%f',channel,y));
  sigma = fsolve(fun,1,optimset('fsolve'));


En = 1;    % set constant energy/channel use
N0 = 2* sigma^2;
Eb = En/r;
SNR = Eb./N0;
SNRdB = 10*log10(SNR);
  