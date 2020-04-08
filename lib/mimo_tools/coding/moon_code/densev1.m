function mulist = densev1(EbN0,Nit,wc,wr,stopmu)
% function mulist = densev1(EbN0,Nit,wc,wr)
%
%Compute examples of density evolution 
%
% EbN0 -- message SNR (not in dB)
% Nit -- number of iterations
% wc -- column weight
% wr -- row weight
% stopmu -- (optional) stop iterations if mu exceeds this value

% Todd K. Moon, May 7, 2004
% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

global psi2
global mlist2

if(length(psi2)==0)  % if not already created
  % Compute and save values of the Psi function
  % (So this only needs to be done once
  fprintf(1,'Please wait: generating Psi function');
  fprintf(1,' (only needs to be done once)\n'); 
  
  mlist = .01:.01:20; 					%only need positive values: psi is odd
  psi = [];
  for m=mlist;
	s2 = 2*m;
    str=sprintf('tanh(x/2).*exp(-(x-(%g)).^2 ./(4*(%g)))/sqrt(4*pi*(%g))',...
	m,m,m);
    F = inline(str);
%    Q = quad(F,m-20*s2,m+20*s2);
    Q = quad(F,m-10*s2,m+10*s2);  % this seems to work as well as 20*s2
	psi = [psi Q];
  end
  mlist2 = [-fliplr(mlist) 0 mlist];
  psi2 = [-fliplr(psi) 0 psi];
end

checkstopmu=0;
if(nargin==5)
  checkstopmu = 1;
end
  

%xlist = 0:.1:8;
%plot(xlist,1-exp(-.4527*(xlist).^.86 + 0.0218),'k'); % plot an appoximation

% To evaluate Psi(x), use interp1(mlist2,psi2,x)
% to evaluate Psi^{-1}(x), use interp1(psi2,mlist2,x)

% Set up the problem parameters
R = 1-wc/wr;

ecs2 = 2*2*R*EbN0; 						% Ec/sigma^2, used in iteration

mu = 0;
mulist = [];
for i=1:Nit
  x = ecs2 + (wc-1)*mu;
%  psi = interp1(mlist2,psi2,x);
  psi = Psi(x,mlist2,psi2);
  psi = psi^(wr-1);
  mu = Psiinv(psi,mlist2,psi2);
  mulist = [mulist mu];
  if(checkstopmu & mu > stopmu) break; end;
end
