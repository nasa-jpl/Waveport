




L = 60;
alpha = pi/3;
beta = pi/4;
gamma = pi/5;
str = [];

tic
[row col Dlmp] = DlmpSparse(L,alpha,beta,gamma,str);
toc

tic
[row2 col2 Dlmp2] = DlmpSparseFast(L,alpha,beta,gamma,str);
toc

c1 = mean(abs(row2-row));
c2 = mean(abs(col2-col));
c3 = mean(abs(Dlmp2-Dlmp));

[c1 c2 c3]