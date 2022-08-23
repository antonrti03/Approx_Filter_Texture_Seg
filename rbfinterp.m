function [f] = rbfinterp(x, options)
tic;
phi       = options.('rbfphi');
rbfconst  = options.('RBFConstant');
nodes     = options.('x');
rbfcoeff  = (options.('rbfcoeff'))';
[dim              n] = size(nodes); %nodes = [cols; rows]
[dimPoints  nPoints] = size(x); %x=[cols; rows]

if (dim~=dimPoints)
  error(sprintf('x should have the same number of rows as an array used to create RBF interpolation'));
end;

f = zeros(1, nPoints);
r = zeros(1, n);

for i=1:1:nPoints
	s=0;
    r =  (x(:,i)*ones(1,n)) - nodes;%?[3;4]*ones(1,10)=[3 3...; 4 4...]
    r = sqrt(sum(r.*r, 1));
%     for j=1:n
%          r(j) =  norm(x(:,i) - nodes(:,j));
%     end
    
     s = rbfcoeff(n+1) + sum(rbfcoeff(1:n).*feval(phi, r, rbfconst));
 
	for k=1:dim
       s=s+rbfcoeff(k+n+1)*x(k,i);     % linear part
	end
	f(i) = s;
end;

if (strcmp(options.('Stats'),'on'))
    fprintf('Interpolation at %d points was computed in %e sec\n', length(f), toc);    
end;
