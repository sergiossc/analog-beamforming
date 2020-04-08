% Plot the Hamming and Gilbert bounds for binary codes

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

% r = n-k
nlist=50:10:1000;
t = 3;
rminlist = [];
rmaxlist = [];
for n = nlist
  rmin = log2(Hammsphere(n,2,t));
  rmax = log2(Hammsphere(n,2,2*t));
  rminlist = [rminlist rmin];
  rmaxlist = [rmaxlist rmax];
end
plot(nlist,rminlist);
hold on;
plot(nlist,rmaxlist,'r');
  