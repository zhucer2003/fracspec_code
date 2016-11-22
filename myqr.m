ccc

n = 20;
m = 1;
seedRNG(0);


% A = rand(n);
% A = tril(triu(A, -m+1), m-1);
% A(1,:) = 1+0*rand(1,n);
% b = rand(n,1);

A = ultraS.diffmat(n,1) + ultraS.convertmat(n, 1, 1);
e = ones(1,n); 
A = [e ; A(1:(end-size(e,1)),:)];
b = rand(n,1).*exp(-(1:n)');
A = full(A);

A0 = A;
b0 = b;
x0 = A0\b0;

R = A;
I = speye(n);
K = n;
bk = ones(n,1);
for k = 1:K-1
    for j = m+1:-1:2
        l = k+j-2;
        
        a1 = R(l+1,k);
        a2 = R(l,k);

        c = a2/sqrt(a1^2+a2^2);
        s = -a1/sqrt(a1^2+a2^2);
        G = I;
        G(l:l+1,l:l+1) = [c -s ; s c];

        bk(k+1) = s*bk(k);
        bk(k) = c*bk(k);

        R = G*R;
        b = G*b;
        R(l+1,k) = 0;
    end
    
    if ( k > 2 )
%         R(k-2,k+m+1:end) = 0;
    end
end

R = R(1:K,1:K);
b = b(1:K);

x1 = R\b;
norm(x1 - x0(1:K), inf)

B = tril(R, 2*m);
p = 0;
x2(K-1:K,1) = B(K-1:K,K-1:K)\b(K-1:K);
for k = K-m-1:-1:1
    s = 1:m+1;
    x2(k) = 1./B(k,k)*(b(k) - B(k,k+s)*x2(k+s) - bk(k)*p);
    p = p + x2(k+m+1)*A0(1,k+m+1);
end
norm(x2- x0(1:K), inf)








