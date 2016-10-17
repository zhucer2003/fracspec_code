function [D, lam, mu] = diffmat(n, m, lam)
if ( nargin == 2 )
    lam = m; m = 1;
end
if ( round(m) ~= m )
    [D, lam, mu] = diffmat05(n, lam);
else
    D = eye(n); mu = lam-.5;
end
for k = 1:floor(m)
    [Dk, lam, mu] = diffmat1(n, lam, mu);
    D = Dk*D;   
end
end

function [D, lam, mu] = diffmat05(n, lam)
    e = (gamma(lam+.5)/gamma(lam)) * ones(n,1);
    D = spdiags([e,e], [0,1], n, n);
    mu = lam - 1;  lam = lam + .5;
end

function [D, lam, mu] = diffmat1(n, lam, mu)
    % d/dx (1+x)^mu*C^{(l)}_n(x)
    e = ones(n, 1);
    if ( mu == 0 )
        D = 2*lam*spdiags(e, 1, n, n);
    else
        v = (mu-lam)./((0:n-1).'+lam);
        D = spdiags(lam*[1+v, 2*e, 1-v], 0:2, n, n);
        mu = mu - 1;
    end
    lam = lam + 1;
end