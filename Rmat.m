function [R, l] = Rmat(n, l)
if ( l == 0 )
    e = ones(n,1);
    R = spdiags(.5*[e,2*e,e], -1:1, n, n);
    R(2,1) = 1;
else
    nn = (0:n-1).';
    e = ones(n,1);
    e1 = .5*(nn+1)./(nn+l);
    e2 = .5*(nn+2*l-1)./(nn+l);
    R = spdiags([e1,e,e2],-1:1,n,n);
end
end