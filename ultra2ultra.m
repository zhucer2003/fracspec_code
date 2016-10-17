function c = ultra2ultra(c, lam1, lam2)

n = length(c) - 1;
c = c./scl(lam1, n);
c = jac2jac(c, lam1-.5, lam1-.5, lam2-.5, lam2-.5);
c = c.*scl(lam2, n);

end

function s = scl(lam, n)
    if ( lam == 0 )
        nn = (0:n-1).';
        s = [1 ; cumprod((nn+.5)./(nn+1))];
    else
        nn = (0:n).';
        s = ( gamma(2*lam) ./ gamma(lam+.5) ) * ...
            exp( gammaln(lam+.5+nn) - gammaln(2*lam+nn) );
    end
end