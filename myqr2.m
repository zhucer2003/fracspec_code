ccc

n = 9;
m = 2;
K = 1;
seedRNG(0);

A = rand(n);
A = tril(triu(A, -m), m);
A(1,:) = rand(1,n);

% A = ultraS.diffmat(n,1) + ultraS.convertmat(n, 1, 1);
% e = ones(1,n); 
% A = [e ; A(1:(end-size(e,1)),:)];

b = rand(n,1).*exp(-(1:n)');
A = full(A);

A0 = A;
b0 = b;
x0 = A0\b0;

R = A;
I = speye(n);
k_max = n;
bk = ones(n,1);
for k = 1:k_max-1
    for j = m+1:-1:2
        l = k+j-2;
        if ( l >= n ), continue, end
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

R = R(1:k_max,1:k_max);
b = b(1:k_max);

x1 = R\b;
norm(x1 - x0(1:k_max), inf)


bk = R(:,end)/A0(1,end);

B = tril(R, 2*m);
p = 0;
x2(k_max-2*m:k_max,1) = B(k_max-2*m:k_max,k_max-2*m:k_max)\b(k_max-2*m:k_max);
x2a = x2;
p = x2(end)*A0(1,end);
for k = k_max-2*m-1:-1:1
    s = 1:2*m;
    x2(k) = 1./B(k,k)*(b(k) - B(k,k+s)*x2(k+s) - bk(k)*p);
    p = p + x2(k+2*m)*A0(1,k+2*m);
end
norm(x2- x0(1:k_max), inf)

% x3 = (B + triu(bk*A0(1,:),2*m+1))\b;
% norm(x3- x0(1:k_max), inf)

[x1 x2]








