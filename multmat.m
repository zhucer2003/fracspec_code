function M = multmat(n, f, lam)

if ( isa(f, 'function_handle') )
    c = chebtech1.vals2coeffs(f(chebpts(n,1))); % T coeffs
    c = ultra2ultra(c, 0, lam);                 % C^(lam) coeffs
else
    c = f(:);
end

% Trim trailing zeros from c (or entries below tolerance):
tol = 100*norm(c, inf)*eps; % TODO: Better tolerance.
idx = find(abs(flipud(c)) > tol, 1, 'first');
c(n-idx+2:end) = [];
c(abs(c) < tol) = 0;
if ( isempty(c) ), M = 0; return, end

% Jacobi matrix:
nn = (0:n)';
J = spdiags(.5*[(nn+1)./(nn+lam), (nn+2*lam-1)./(nn+lam)], [-1,1], n, n);

% Initialise recurrence:
Cm1 = sparse(0); C = speye(n);
M = c(1)*C;

% Recurrence relation:
for n = 0:length(c)-2
    Cp1 = (2*(n+lam)/(n+1))*(J*C) - ((n+2*lam-1)/(n+1))*Cm1;
    Cm1 = C; C = Cp1;
    if ( c(n+2) ~= 0 )
        M = M + c(n+2)*C;
    end
end

end