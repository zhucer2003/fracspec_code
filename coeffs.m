function c = coeffs(f, n, lam)
if ( isnumeric(f) )
    c = zeros(n,1); c(1) = f;
else
    c = ultracoeffs(chebfun(f), n, lam);
end
end