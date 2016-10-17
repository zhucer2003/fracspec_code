function [Q, lam, mu] = cumsummat(n, m, lam)
    if ( nargin == 2 )
        lam = m; m = 1;
    end
    Q = 1;
    Q1 = cumsummat1(n, lam);
    for k = 1:m
        Q = Q1*Q;
    end
    if ( round(m) ~= m )
        Q = cumsummat05(n, lam)*Q;
        lam = lam-floor(lam)+.5;
    end
    mu = lam - .5;
end

function [Q, lam, mu] = cumsummat05(n, lam)
    if ( lam == .5 )
        v = (2/sqrt(pi)) ./ (2*(0:n-1).'+1);
        Q = spdiags([v,-v], [0,1], n, n);
        mu = .5; lam = 1;
    elseif ( lam == 1 )
        v = (sqrt(pi)/2) * ones(n,1);
        Q = spdiags([v,v], [-1,0], n, n);
        mu = lam;  lam = lam + 0.5;
    end
end

function [Q, lam, mu] = cumsummat1(n, lam)
    if ( lam == .5 )
        v = 1 ./ (2*(0:n-1).'+1);
        Q = spdiags([v,-v], [-1,1], n, n); Q(1,1) = 1;
        mu = 0;
    elseif ( lam == 1 )
        nn = (0:n)';
        v = 1 ./ (2*nn+1);  v2 = 2 ./ (4*nn.^2-1);
        Q = spdiags([v(2:end),v2(2:end),-v(1:end-1)], [-1:1], n, n);
        mu = 0.5;
    end
end