% Plot the probability of error for a Hamming code as a function of SNR

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only


SNR = 0:1:12;  % Eb/N0, in dB
EbN0 = 10.^(SNR/10);

Pnocode = qf(sqrt(2*EbN0));
clf;

% Stuff for Hamming code
% for (7,) code:
% A(x) = 1/(n+1)( (1+x)^n + n(1-x)(1-x^2)^(n-1)/2)
n = 7;
k = 4;
EcN0 = EbN0*k/n;
Pen = qf(sqrt(2*EcN0));

% Build the weight enumerator
xp1 = [1 1];   % 1+x
xm1 = [-1 1];  % 1-x;
xm2 = [-1 0 1]; % 1-x^2
A1 = xp1;
for i=1:n-1    % compute (1+x)^n
  A1 = conv(A1,xp1);
end
A2 = xm2;
for i=1:(n-1)/2-1  % compute (1-x^2)^(n-1/2)
  A2 = conv(A2,xm2);
end
A = (A1 + n*conv(xm1,A2))/(n+1);
dmin = 3;  % for Hamming codes
t = (dmin-1)/2;

j1 = 0;
for p = Pen
  j1 = j1+1;

  % compute the P_k^j functions of (10-11) of Wicker
  for j=dmin:n
	for k1=0:t
	  P(k1+1,j) = 0;
	  for r=0:k1
		P(k1+1,j) = P(k1+1,j) + nchoosektest(j,k1-r)*nchoosektest(n-j,r)* ...
			p^(j-k1+2*r)*(1-p)^(n+j+k1-2*r);
	  end
	end
  end
  % Compute P(E), the probability of decoder error:
  PE(j1) = 0;
  for j=dmin:n
	s = 0;
	for k1=0:t
	  s = s + P(k1+1,j);
	end
	PE(j1) = PE(j1) + A(j+1)*s;
  end
end



% Compute B_j
Par = [1 1 0 1
	   1 0 1 1
	   0 1 1 1];
H = [eye(n-k) Par];
G = [Par' eye(k)];
B = zeros(n,1);
for i=1:2^k-1
  m = dec2bin(i,k)-'0';  % get message vector
  wm = sum(m);
  c = mod(m*G,2);
  wc = sum(c);
  B(wc) = B(wc) + wm;
end

j1 = 0;
for p = Pen
  j1 = j1+1;

  
  % compute the P_k^j functions of (10-11) of Wicker
  for j=dmin:n
	for k1=0:t
	  P(k1+1,j) = 0;
	  for r=0:k1
		P(k1+1,j) = P(k1+1,j) + nchoosektest(j,k1-r)*nchoosektest(n-j,r)* ...
			p^(j-k1+2*r)*(1-p)^(n+j+k1-2*r);
	  end
	end
  end

  % Compute Pb, the probability of decoder error:
  Pb(j1) = 0;
  for j=dmin:n
	s = 0;
	for k1=0:t
	  s = s + P(k1+1,j);
	end
	Pb(j1) = Pb(j1) + B(j)*s;
  end
  Pb(j1) = Pb(j1)/k;
end

semilogy(SNR,Pnocode);
hold on
semilogy(SNR,Pen,'r--');
%semilogy(SNR,PE,'r:');
%semilogy(SNR,PE/k,'g:');
semilogy(SNR,Pb,'-.');
legend('BER with no coding','BER with coding','Decoded BER')

cmpdb = 10;
idx = find(SNR==cmpdb);
Pb(idx);
Pcoded = Pb(idx);

semilogy([cmpdb,cmpdb],[Pnocode(idx), Pb(idx)],':k');

SNReff = interp1(log(Pnocode),SNR,log(Pcoded),'linear');
fprintf(1,'Uncoded: %g  coded: %g  SNReff=%g  gain:%g\n',Pnocode(idx),...
	Pb(idx),SNReff,SNReff-cmpdb);
codegain = SNReff - cmpdb;
semilogy([cmpdb,SNReff],[Pb(idx),Pb(idx)],'k');
semilogy([cmpdb,cmpdb],[Pb(idx)/2,Pb(idx)*2],'k');
semilogy([SNReff,SNReff],[Pb(idx)/2,Pb(idx)*2],'k');

text((SNReff+cmpdb)/2,Pb(idx)/3,'Coding gain','HorizontalAlignment','center');
xlabel('SNR(dB)');
ylabel('Probability of bit error');
print -dps hamcode74pe.ps
