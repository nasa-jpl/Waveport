



%% cubic

M = 1;
N = 1;

a = 4*rand(M,N);
b = rand(M,N);
c = -10*rand(M,N);
d = rand(M,N);

coef1 = [a b c d];
[x1 x2 x3] = cubicroots(a,b,c,d)

[[x1 x2 x3].' roots(coef1)]


%% quartic


M = 1;
N = 1;

rng(2);
a = 4*rand(M,N);
b = 0*rand(M,N);
c = -10*rand(M,N);
d = 0*rand(M,N);
e = 0*rand(M,N);

coef1 = [a b c d e];

[x1 x2 x3 x4] = quarticroots(a,b,c,d,e);

for m=1:M,
for n=1:N,
    roots([a(m,n) b(m,n) c(m,n) d(m,n) e(m,n)])
    pause
end
end


[[x1 x2 x3 x4].' roots(coef1)]













