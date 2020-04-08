% plot and compare the capacity of the AWGNC and the BAWGNC

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

snr = 0.001:0.25:6.001;  % signal to noise ratio of Ec/sigma^2
cawgnc = [];
cbawgnc = [];
cbsc = [];
for s = snr
  cawgnc = [cawgnc cawgnc2(s)];
  cbawgnc = [cbawgnc cbawgnc2(s)];
  p = qf(sqrt(s));
  cbsc = [cbsc 1-h2(p)];
end



clf
plot(snr,cawgnc);
hold on;
plot(snr,cbawgnc,'--');
plot(snr,cbsc,'-.');
xh = xlabel('E_c/\sigma^2');
yh = ylabel('Capacity (bits/channel use)');
axis([0 6 0 1.5])
drawnow
lh = legend('AWGNC Capacity', 'BAWGNC Capacity','BSC capacity');
set(gca,'fontsize',15);
set(xh,'fontsize',15);
set(yh,'fontsize',15);
set(lh,'fontsize',15);
grid on
input('press return to save');
print -dps capcmp.ps
