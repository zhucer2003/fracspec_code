ccc
mu = .5;
Nmax = 50;

trueSoln = @(x) exp(1+x).*erfc(sqrt(1+x))-1;

for N = [Nmax+21 2:Nmax]
    
    [D, x] = cheb_mu(N, mu);
    I = eye(N);

    A = I + D;
   
    A = A(2:end,2:end);
    x = x(2:end);

    v = A\(-1+0*x);

    err(N) = norm(v(end) - trueSoln(x(end)), inf);

    figure(1)
    plot([-1 ; x], [0 ; v]), shg

end

figure(2)
N = 1:N;
semilogy(N, err(N), '-', N, .01*exp(-sqrt(N)), '-');

figure(3)
loglog(N, err(N), '-', N, .01*exp(-sqrt(N)), '-');

alignfigs

