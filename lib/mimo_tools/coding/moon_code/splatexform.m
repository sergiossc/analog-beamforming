function latexform(fid,X,nohfill)
% function latexform(fid,X,[nohfill])
% display a matrix X in latex form
[n,m] = size(X);
if(nargin == 3)
  nohfill = 1;
else
  nohfill = 0;
end
fprintf(fid,'\\begin{bmatrix}');
for i=1:n
  for j=1:m
	if(X(i,j))
	  if(nohfill)
		fprintf(fid,'%.2g ',X(i,j));
	  else
		fprintf(fid,'\\hfill%.2g ',X(i,j));
	  end
	end
    if(j ~= m) fprintf(fid,'& ');
    end
  end
  if(i ~= n)
    fprintf(fid,'\\\\\n');
  end
end
fprintf(fid,'\\end{bmatrix}\n');