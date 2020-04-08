% Todd K. Moon, September 30, 2003

% Plot the capacity for a Gaussian and Binary AWGN channel

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

figure(1);
clf;

markerlist = 'ox+*sd';
channellist{1} = 'cbawgnc';
channellist{2} = 'cawgnc'; 		% list of channels to use
channelname = {'CBAWGNC','AWGNC'};
leglist = {};
linestyle={'-','--'};

Rlist = [1/4 1/2 3/4 .9]; 				% list of rate to plot for
Pb = logspace(-2,-6,20);  % probabilities of error betwen 10^2 and 10^-6
c1 = 0;
d1 = 0;
llist = [];

for c1 = 1:length(channellist) 				% different channels to use
  channel = channellist{c1};
  r1 = 0;
  for r = Rlist
	r
	r1 = r1+1;
	d1 = d1+1;
	leglist{d1} = sprintf('%s R=%.4g',channelname{c1},Rlist(r1));
	sigmalist = zeros(1,length(Pb));
	i = 0;
	for p = Pb
	  i = i+1;
	  y = r*(1-h2(p)); 	% find rate after compression with distortion p
	  fun = inline(sprintf('%s(sigma)-%f',channel,y));
	  sigmalist(i) = fsolve(fun,1,optimset('fsolve'));
	  % solve for SNR which is equal to y
	end
	
	Ec = 1; 							% set constant energy/channel use
	N0 = 2* sigmalist.^2;  % convert variance to N0
	Eb = Ec/r;
	SNR = Eb./N0;          % convert to Eb/N0
	SNRdB = 10*log10(SNR); % convert to dB
	llist = [llist semilogy(SNRdB,Pb,[linestyle{c1} markerlist(r1)])]; 
	hold on;
	drawnow;
  end
end
lh = legend(llist,leglist);
xh = xlabel('E_b/N_0 (dB)');
yh = ylabel('P_b');
set(gca,'fontsize',15);
set(xh,'fontsize',15);
set(yh,'fontsize',15);
set(lh,'fontsize',15);
input('press return to print');
print -dps plotcbawgn2.ps
