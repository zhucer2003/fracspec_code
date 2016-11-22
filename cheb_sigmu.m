function [D, x] = cheb_sigmu(N, sigma, mu)

% N = 15;
% clc
% mu = .125;
% sigma = .5;

% Chebyshev grid:
x = -cos((0:(N-1))*pi/(N-1))';

%% Compute beta (v3)
V = zeros(N);
V(:,1) = 1;    
V(:,2) = x - mu;   
for k = 2:N-1 % Jacobi recurrence.
    km1 = k - 1;
    V(:,k+1) = (2*(k - .5)*x.*V(:,k) - (km1 - mu^2/km1)*V(:,k-1)) / k;
end
bet = inv(V);

cl = ceil(sigma-mu);

foo = @(n,k) gamma(n+1)/(gamma(k+1)*gamma(n-k+1));
for n = 1:N
    for q = 1:n-1
        b(n,q) = (-1)^(n+q-1)*(.5)^q*foo(n-1+q,q)*foo(n-1+mu,n-1-q)*gamma(q+mu+1)/gamma(q+mu-sigma+1);
    end
end

D = zeros(N);
for i = 2:N
    for j = 2:N
        for n = 1:N
            tmp = 0;
            for q = max(cl,1):(n-1)
                tmp = tmp + b(n,q)*(x(i)+1)^(q+mu-sigma);
            end
            
            D(i,j) = D(i,j) + bet(n,j)*tmp;
        end
        D(i,j) = D(i,j) / (1+x(j))^mu ;
    end
end


%%
% Truncate boundary point:
% D = D(2:N,2:N);
% x = x(2:end);

end