clc

a = .5;

f = chebfun(@(x) x.^(a-1).*ml(-x.^a, a, a), [0, 1], 'exponents', [-.5, 0], 10000)
f = chebfun(@(x) (1+x).^(a-1).*ml(-(1+x).^a, a, a), 'exponents', [-.5, 0], 10000)

figure(1)
plot(f)
figure(2)
plot(diff(f, a) + f), shg
alignfigs