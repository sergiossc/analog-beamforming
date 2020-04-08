% Make plots of density evolution stuff

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

if 0
figure(1);
clf;

mu172 = densev1(10^(1.72/10),800,4,6,1);
plot(mu172)
drawnow;
% mu = .3155

hold on;

% mu174 = densev1(10^(1.74/10),800,4,6,1);
% plot(mu174)
% drawnow;

mu175 = densev1(10^(1.75/10),800,4,6,1);
plot(mu175)
drawnow;
% mu=0.3564

mu176 = densev1(10^(1.76/10),800,4,6,1);
plot(mu176)
drawnow;
% mu=0.3832

% mu1763 = densev1(10^(1.763/10),800,4,6,1);
% plot(mu1763)
% drawnow;
% % doesn't break away

mu1764 = densev1(10^(1.764/10),800,4,6,1);
plot(mu1764)
drawnow;
% breaks away at about 500

mu1765 = densev1(10^(1.765/10),800,4,6,1);
plot(mu1765)
drawnow;
% breaks away at about 265

mu177 = densev1(10^(1.77/10),800,4,6,1);
plot(mu177)
drawnow;
% breaks away at about 120

mu178 = densev1(10^(1.78/10),800,4,6,1);
plot(mu178)
drawnow;
% breaks away at about 72

mu180 = densev1(10^(1.80/10),800,4,6,1);
plot(mu180)
drawnow;
% breaks away at about 45

axis([0 800 0 1]);
hx = xlabel('iteration l');
hy = ylabel('\mu^{[ l ]}');
set(gca,'fontsize',15);
set(hx,'fontsize',15);
set(hy,'fontsize',15);
print -deps densev1.eps

end


figure(2); 
clf

[x172,pg172] = plotgauss(mu172(end),2*mu172(end));
plot(x172,pg172,'--');
drawnow

hold on;

mu1764 = densev1(10^(1.764/10),800,4,6,80);
for m=507:length(mu1764)
  [x1764,pg1764] = plotgauss(mu1764(m),2*mu1764(m));
  plot(x1764,pg1764)
  drawnow
end
axis([-10 100 0 0.55])
xh = xlabel('\lambda^{[ l ]}');
yh = ylabel('p(\lambda^{[ l ]})');
set(gca,'fontsize',15);
set(xh,'fontsize',15);
set(yh,'fontsize',15);
print -deps densev2.eps
