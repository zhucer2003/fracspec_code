function [u, A, sol, rhsOut] = example10(n)

bc = [1 ; 0];
e = 0;
f = 0;
sol = @(x) NaN*x;
    
% Initialise block operators
init
idx = [1 ; 2 ; idx+2];

% Construct the operator for this example:
A1 = [e1;z];
A2 = QQ(.5)*[S(.5)\(D(1,.5)*e2) ; z]+[e2;z];
A3 = II + QQ(1.5) + QQ(2);
A = [A1, A2, A3];
% Construct the rhs:
rhs = [coeffs(e, n, .5) ; coeffs(f, n, 1)];

% Construct and append boundary conditions:
ee = ones(1,n); ee2 = ee; ee2(2:2:end) = -1; ee3 = 1:n;
BC = [ee2, z' ; ee, sqrt(2)*ee3]*QQ(2);
A = [[1 -1 ; 1 1] BC ; A];
rhs = [bc ; rhs];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Solve for v
v(idx,1) = A\rhs;
% v(idx,1) = mysolve(A, rhs, 2);

% Reconstruct u:
u(1:2*n,1) = QQ(2)*v(3:end);
u(1:2) = u(1:2) + v(1:2);

rhsOut = @(x) e(x) + sqrt(1+x).*f(x);

end