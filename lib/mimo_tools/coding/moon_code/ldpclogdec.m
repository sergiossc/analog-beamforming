function x = loggaldecode(A,r,Nloop,Lc) 
% function x = loggaldecode(A,r,Nloop,Lc) 
%
% Do log-likelihood decoding on a low-density parity check code
% A = parity check matrix
% r = received signal vector
% Nloop = number of iterations
% Lc = channel reliability

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

[M,N] = size(A);
clear Nl Ml
for m=1:M Nl{m} = []; end;
for n=1:N Ml{n} = []; end;
% Build the sparse representation of A using the M and N sets
for m=1:M
  for n=1:N
	if(A(m,n))
	  Nl{m} = [Nl{m} n];
	  Ml{n} = [Ml{n} m];
	end
  end
end
idx = find(A ~= 0);  % identify the "sparse" locations of A
% The index vector idx is used to emulate the sparse operations
% of the decoding algorithm.  In practice, true sparse operations
% and spare storage should be used.

% Initialize the probabilities
eta = zeros(M,N);
lasteta = zeros(M,N);
lambda = Lc*r;
% fprintf(1,'lambda[0]:'); splatexform(1,lambda,1);

for loop = 1:Nloop
  fprintf(1,'loop=%d\n',loop);
  for m = 1:M 			  % for each row (check)
	for n=Nl{m} % work across the columns ("horizontally")
	  pr = 1;
	  for np= Nl{m}
		if(np == n) continue; end;
		pr = pr*tanh((-lambda(np)+lasteta(m,np))/2); % accumulate the product
	  end
	  eta(m,n) = -2*atanh(pr);
	end
  end
  lasteta = eta;   % save to subtract to obtain extrinsic for next time around
  % fprintf(1,'eta:'); splatexform(1,eta,1);

  for n=1:N 			        % for each column (bit)
	lambda(n) = Lc*r(n);
	for m = Ml{n}	% work down the rows ("vertically")
	  lambda(n) = lambda(n) + eta(m,n);
	end
  end
%   fprintf(1,'lambda:'); splatexform(1,lambda,1);
%   p = exp(lambda) ./ (1+exp(lambda));  % needed only for comparison purposes!
%   fprintf(1,'p:'); splatexform(1,p,1);

  x = lambda >= 0;  % compute decoded bits for stopping criterion
  z1 = mod(A*x',2)';

%  fprintf(1,'x: ');latexform(1,x,1);
%  fprintf(1,'z: ');latexform(1,z1,1);
  if(all(z1==0)) break; end
end  % end for loop

if(~all(z1==0))
  fprintf(1,'Decoding failure after %d iterations',Nloop);
end