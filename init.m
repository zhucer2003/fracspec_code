
% Basic:
Z = sparse(n,n);
I = speye(n);
II = [I Z ; Z I];
idx = [1:n ; n+1:2*n]; idx = idx(:);

e1 = sparse([1 ; zeros(n-1,1)]);
e2 = sparse([0 ; 1 ; zeros(n-2,1)]);
z = sparse(zeros(n,1));

% Conversion:
S = @(lam) Smat(n, lam);
R = @(lam) Rmat(n, lam);
EE5 = @(m) [S(m+.5) Z ; Z R(m+1)];
EE1 = @(m) [I Z ; Z S(m)];
EE = @(m) (round(m)==m)*EE1(m) + (round(m)~=m)*EE5(m-.5);

% Jacobi matrix:
J = @(lam) spdiags(.5*[((0:n)+1)./((0:n)+lam);  ((0:n)+2*lam-1)./((0:n)+lam)]', [-1,1], n, n);

% Integration:
Q05 = @(lam) cumsummat(n, .5, lam);
QQ05 = [Z Q05(1) ; Q05(.5) Z];
QQ = @(m) QQ05^(2*m);

% Differentiation:
D = @(m,lam) diffmat(n, m, lam);
DD1 = @(m) [D(m,.5) Z ; Z D(m,1)];
DD5 = @(m) [Z D(m,1) ; D(m,.5) Z];
DD = @(m) (round(m)==m)*DD1(m) + (round(m)~=m)*DD5(m);

% Connection:
C = @(a,b) triu(ultra2ultra(eye(n),a,b));

% Multiplication:
X = blkdiag(J(.5), J(1));
Pi1 = @(f, lam) multmat(n, f, lam);
% Pi = @(lam, r, s) [Pi1(r,lam), C(lam+.5,lam)*Pi1(@(x)(1+x).*s(x), lam+.5) ;
%                    Pi1(s,lam+.5)*C(lam,lam+.5), Pi1(r, lam+.5)];
Pi = @(lam,r) blkdiag(Pi1(r, ceil(lam) + .5), Pi1(r, floor(lam) + 1));
                    