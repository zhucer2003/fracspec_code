function [u, A] = example11(n, ep)

% Initialise block operators
init
% Construct the operator for this example:
% A = ep*1i^1.5*DD(1.5) - Pi(1.5,@(x) x)*EE(1.5)*EE(1)*EE(.5);
% A = ep*1i^1.5*DD(1.5) - EE(1.5)*EE(1)*EE(.5)*X;
A = ep*1i^1.5*DD5(1.5) - EE5(1)*EE1(1)*EE5(0)*X;

% Construct the rhs:
rhs = zeros(2*n, 1);

% Construct and append boundary conditions:
BCR = [ones(1,n), sqrt(2)*(1:n)]; % Right Dirichlet
BCL = [ones(1,n), zeros(1,n)]; BCL(2:2:n) = -1; % Left Dirichlet

%%
bc = [0 ; 1];
BC = [BCL ; BCR]; 

nbc = length(bc);

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);
BC = BC(:,idx);

A = [BC ; A(1:end-nbc,:)];
rhs = [bc ; rhs(1:end-nbc)];

% Solve:
u = [];
u(idx,1) = mysolve(A, rhs, 2);

end
