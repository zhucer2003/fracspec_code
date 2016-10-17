function [u, A, sol] = example6(n)

bc = 1;
e = 0;
f = 0;
sol = @(x) NaN*x;
    
% Initialise block operators
init
idx2 = idx(1:end-1);

% Construct the operator for this example:
A = DD(1) + EE(1)*(DD(0.5) + EE(.5));
% Construct the rhs:
rhs = [coeffs(e, n, .5) ; coeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Construct and append boundary conditions:
BC = [ones(1,n), sqrt(2)*(1:n)]; % Right Dirichlet
BC(2:2:end) = -BC(2:2:end); BC(n+1:end) = 0; % Left Dirichlet
A = [BC(idx) ; A];
rhs = [bc ; rhs];

% Solve:
u(idx,1) = A\rhs;

end