% a Generalized distributive law message passing algorithm
% This actually deals with junction trees instead of factor
% graphs, so it is somewhat more general.

% See: Aji and McEliece, The Generalized Distributive Law,
% IEEE Trans. Info. Th. v. 46, Mar. 2000, pp. 325--343

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

% input:
% 1) A connection list such as indicating connections between nodes
%    e.g., C{1} = [1, 11 12 15] ...  C{11} =  [11, 1,2,3,6,7,10]
% 2) a list of function names, one for each node
%    in the cell array F
%    e.g. F{1} = 'fxy'  F{2} = 'fxy', etc.
% 3) a list of auxilliary function arguments in the cell array Y
%    e.g. Y{1} = y(1)... Y{11} = [], etc.
% 4) A list of variables at each node in S
%    e.g. S{1} = [1]  S{2} = [2] ... S{11} = [1 2 3 6 7 10]
% 5) A number of iterations Niter
% 6) The number of variables Nvar
% 7) The number of distinct values of each variable in Q
% 8) The number of nodes for which to compute/print the beliefs in Nin
% 
% For an example, see testgdl2.m

Nnode = length(C); 						% number of nodes
clear mu;
% Initialize all the messages to 1
for n = 1:Nnode
  i = C{n}(1);   % from node
  for j = C{n}(2:end)   % to node
	 SiINTSj = intersect(S{i},S{j});
	 % Evaluate for each variable in SiINTSj
	 Neval = length(SiINTSj);
	 xmax = Q^Neval-1;
	 for xsiintsjctr = 0:xmax  % for each value in the table
	   mu{i}{j}(xsiintsjctr+1) = 1;
	 end
  end
end

% cycle1

dilist = [];

for Ni=1:Niter
  % Do the message passing
  for n = 1:Nnode
	i = C{n}(1);   % from node
	for j = C{n}(2:end) 				% to node
	   SiINTSj = intersect(S{i},S{j});
	   SiDIFFSj = setdiff(S{i},S{j});
	   % Evaluate for each variable in SiINTSj
	   Neval1 = length(S{i});
	   Neval2 = length(SiINTSj);
	   xmax = Q^Neval2;
	   Nsum = length(SiDIFFSj);
	   xsummax = Q^Nsum;
	   
	   for xsiintsjctr = 0:xmax-1  % for each value in the table
		 xsiintsj = numtoQar(xsiintsjctr,Q,Neval2,SiINTSj,Nvar);
		 mu{i}{j}(xsiintsjctr+1) = 0;

	     % sum over all variables in SiDIFFSj
		 for xsidifsjctr = 0:xsummax-1
		   xsidifsj = numtoQar(xsidifsjctr,Q,Nsum,SiDIFFSj,Nvar);
		   xsi = xsidifsj + xsiintsj;
		   alphai = feval(F{i},xsi(S{i}),Y{i});
	       if(alphai)  % for many problems, alphai is frequently 0
			 for k = C{n}(2:end) 		% vk adj vi
			   if(k == j) continue; end;
			   SkINTSi = intersect(S{i},S{k});
			   xskintsi = numtonum(xsi,SkINTSi,Q);
			   alphai = alphai*mu{k}{i}(xskintsi+1);
		     end
		   end
		   mu{i}{j}(xsiintsjctr+1) = mu{i}{j}(xsiintsjctr+1) + alphai;
		 end  % end sum over all variables
	   end % end for each value in table
	   mu{i}{j} = mu{i}{j}/sum(mu{i}{j});
% fprintf(1,'mu{%d}{%d}=',i,j);  fprintf(1,'%g ',mu{i}{j}); fprintf(1,'\n');
	end % end for j
  end % end for i

  if(Ni > 1)  % if not the first time around
	sigmaold = sigma;
  end

% if(Ni > 0)
% cycle1
% end

  % Now compute the beliefs
  for n = 1:Nin
	i = C{n}(1);   % from node
	Neval1 = length(S{i});
    xmax = Q^Neval1;
	for xsictr = 0:xmax-1
	  xsi = numtoQar(xsictr,Q,Neval1,S{i},Nvar);
	  alphai = feval(F{i},xsi(S{i}),Y{i});
	  if(alphai) % for many problems, alphai is frequently 0
		for k = C{n}(2:end) 			% vk adj vi
		  SkINTSi = intersect(S{i},S{k});
		  xsiintsi = numtonum(xsi,SkINTSi,Q);
		  alphai = alphai*mu{k}{i}(xsiintsi+1);
		end
	  end
	  sigma{i}(xsictr+1) = alphai;
	end  % for xsictr
	sigma{i} = sigma{i}/sum(sigma{i});  % normalize to be a belief
fprintf(1,'sigma{%d}=',i);  fprintf(1,'%g ',sigma{i}); fprintf(1,'\n');
  end % for i
  
  % See how the beliefs have changed
  if(Ni > 1)
	di = 0;
	for n=1:Nin
	  di = di + norm(sigma{n} - sigmaold{n});
	end
	fprintf(1,'di=%g\n',di);
	dilist = [dilist di];
  end


  if(Ni < Niter)
    x = input('Press return');
  end


end % end for Ni
		  
