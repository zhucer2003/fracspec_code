function [S, l] = Smat(n, l)
% C^{(l)}_n(x) --> C^{(l+1)}_n(x)
if ( l == 0 )
    e = ones(n,1);
    S = spdiags(.5*[e,-e], [0,2], n, n);
    S(1) = 1;
else
    e = l./((0:n-1)'+l);
    S = spdiags([e,-e], [0,2], n, n);
end
l = l + 1;
end