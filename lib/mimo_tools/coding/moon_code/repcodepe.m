% Plot the probability of error for repetition codes as a function of SNR

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

SNR = 0:1:12;  % Eb/N0, in dB
EbN0 = 10.^(SNR/10);
Pnocode = qf(sqrt(2*EbN0));


n = 3;         % repetition code length
EcN0 = EbN0/n;

clf;
semilogy(SNR,Pnocode);
hold on;
Pen = qf(sqrt(2*EcN0));
semilogy(SNR,Pen,'--o');
Pdec = 0;
for k=(n+1)/2:n
  Pdec = Pdec + nchoosek(n,k)*Pen.^k .* (1-Pen).^(n-k);
end
semilogy(SNR,Pdec,':o');

n = 11;
EcN0 = EbN0/n;

Pen = qf(sqrt(2*EcN0));
semilogy(SNR,Pen,'--v');
for k=(n+1)/2:n
  Pdec = Pdec + nchoosek(n,k)*Pen.^k .* (1-Pen).^(n-k);
end
semilogy(SNR,Pdec,':v');

legend('BER with no coding','BER with coding, n=3','Decoded BER, n=3',...
	   'BER with coding, n=11', 'Decoded BER, n=11')

xlabel('SNR(dB)');
ylabel('Probability of bit error');
print -dps repcodeps.ps