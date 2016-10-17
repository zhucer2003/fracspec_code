function x = mysolve(A, b, m)
% Fast solution of A*x = b where A is  banded + m dense rows via Schur
% complement factorisation.

n = size(A);

i1 = 1:m;
i2 = m+1:n;
i3 = 2:m+1;

c = A(i2,i2)\[b(i2), A(i2,i1)];
x = (A(i1,i1) - A(i1,i2)*c(:,i3)) \ (b(i1) - A(i1,i2)*c(:,1));
x = [x ; c*[1; -x]];

% norm(M\rhs - x, inf)

end