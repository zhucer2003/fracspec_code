function c = coeffs(f, n, lam)
if ( isnumeric(f) )
    c = zeros(n,1); c(1) = f;
else
    F = chebfun(f);
    F = simplify(F);
    c = ultracoeffs(F, n, lam);
end
end