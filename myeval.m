function uu = myeval(u, x, flag)
n = length(u)/2;
if ( nargin < 3 ), flag = 0; end
if ( flag ) 
    digits = 40;
    u = vpa(u, digits);
end
uu = clenshawP(x, u(1:n), flag) + sqrt(1+x).*clenshawU(x, u(n+1:2*n));
end

function y = clenshawU(x, c)
% Clenshaw scheme for scalar-valued functions.
bk1 = 0*x;
bk2 = bk1;
n = size(c,1)-1;
for k = n:-1:1
    bk = c(k+1) + 2*x.*bk1 - bk2;
    bk2 = bk1;
    bk1 = bk;
end
y = c(1) + 2*x.*bk1 - bk2;
end

function y = clenshawP(x, c, flag)
% Clenshaw scheme for scalar-valued functions.
digits = 30;
bk1 = 0*x; 
bk2 = bk1;
n = size(c,1)-1;
for k = n:-1:1
    if ( flag )
        k_ = vpa(k, digits);
    else
        k_ = k;
    end
    bk = c(k+1) + (2*k_+1)/(k_+1)*x.*bk1 - (k_+1)/(k_+2)*bk2;
    bk2 = bk1;
    bk1 = bk;
end
y = c(1) + x.*bk1 - .5*bk2;
end










