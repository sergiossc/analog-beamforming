% Convert the sigma of [Chung,Richardson, Urbanke, 2001] into Eb/N0

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

sigma = [.8747
	     .8323
		 .791
		 1.0003
		 1.0035
		 1.2517
		 .7440
		 .7051
		 .6297];
R = [.5;
	 .5;
	 .5;
	 .4;
	 1/3;
	 .25;
	 .6;
	 2/3;
	 .75];

EbN0 = 1 ./ (2*R.* (sigma.^2));
EbN0dB = 10*log10(EbN0);
