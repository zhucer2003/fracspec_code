close all

set(0, 'DefaultLineLinewidth', 2)

nn = ceil(logspace(2,2.61,100));
nn = [nn(1) nn];
err2 = 0*nn; err1 = err2;

xstar = linspace(-1,1,101);
xstar = .5;

ep = .0001;


% Error 2:
nMax = nn(end);
n11 = ceil(1.1*nMax);
u2 = example10(n11, ep);
sol = myeval(u2, xstar);

loopnum = 1;
t = zeros(nMax, 1);

for n = nn

    n
    % Solve for each n:
    tic
    for loop = 1:loopnum
        [u, A] = example10(n, ep);
        tmp = toc;
    end
    t(n) = toc/loopnum;
    
    % Error 1:
    err1(n) = norm(myeval(u, xstar) - sol, inf);
end

%%
% Plotting:

figure(1) % Solution
xx = linspace(-1, 1, 10001);
uu = myeval(u, xx);
plot(xx, real(uu), xx, imag(uu));
% ylim([0, 1.1]), grid on
drawnow, shg, pause(eps)
legend('real', 'imag')
print -depsc2 ../fracspec/paper/figures/example10a

%%
figure(5)
N = chebop(@(x,u) ep/10/sqrt(2)*diff(u,2) + x.*u, [-1 1]); 
N.lbc = 0; N.rbc = 1; 
plot(N\0), shg

print -depsc2 ../fracspec/paper/figures/example10a_zoom
alignfigs


%%
figure(2) % Error
% semilogy(nn, err1(nn), '-', nn, err2(nn), '-');
semilogy(nn, err1(nn));
xlim([0, n])
% ylim([1e-16, 1e1])
grid on
drawnow, shg, pause(eps)
% legend('pointwise', 'coefficients')
title('Fractional Airy, eps = 0.0001')
print -depsc2 ../fracspec/paper/figures/example10c

figure(3) % Spy
spy(A)
drawnow, shg, pause(eps)
print -depsc2 ../fracspec/paper/figures/example10d

%%
figure(4)
loglog(nn, t(nn),'-');
hold on
nn2 = nn(ceil(length(nn)/2)+1:end);
plot(nn, 1.1*nn./nn(end)*t(nn(end)), ':', 'color', [1 1 1]*.6);
hold off
drawnow, shg, pause(eps)
ylim([10^(-2.5), 1e0])
grid on
print -depsc2 ../fracspec/paper/figures/example10e

%%
alignfigs
