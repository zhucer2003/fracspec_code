function [u, A, sol, rhsOut] = example9(n)

bc = 1;
e = @(x) 0*x; f = @(x) 0*x;
sol = @(x) exp(1+x).*erfc(sqrt(1+x));
    
% Initialise block operators
init
idx = [1 ; idx+1];
    
% Construct the operator for this example:
A1 = [e1;z];
A2 = QQ(1) + QQ(.5);
A = [A1, A2];
% Construct the rhs:
rhs = [coeffs(e, n, .5) ; coeffs(f, n, 1)];

% Construct and append boundary conditions:
ee = ones(1,n); ee2 = ee; ee2(2:2:end) = -1;
BC = [ee2, z']*QQ(1);

A = [1 BC ; A];
rhs = [bc ; rhs];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);
    
% Solve for v
% v(idx,1) = A\rhs;
v(idx,1) = mysolve(A, rhs, 1);

% Reconstruct u:
u(1:2*n,1) = QQ(1)*v(2:end);
u(1) = u(1) + v(1);

rhsOut = @(x) e(x) + sqrt(1+x).*f(x);

end