function [D, x] = cheb_mu(N, mu)
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

%% Compute P(x_k) (v2)
P = zeros(N);
P(:,1) = 1;    
P(:,2) = x;   
for k = 2:N-1 % Legendre recurrence.
    P(:,k+1) = (2-1/k)*x.*P(:,k) - (1-1/k)*P(:,k-1);
end

%% Compute scl (avoid overflow in gamma function)
% NN = (1:N); % scl = gamma(NN.'+mu)./gamma(NN.')
scl = zeros(N,1);
scl(1) = gamma(1+mu);
for n = 2:N
    scl(n,1) = (1+mu/(n-1))*scl(n-1,1);
end

%%
D = P * spdiags(scl, 0, N, N) * bet * spdiags(1./(1+x).^mu, 0, N, N);

%%
% Truncate boundary point:
% D = D(2:N,2:N);
% x = x(2:end);

end