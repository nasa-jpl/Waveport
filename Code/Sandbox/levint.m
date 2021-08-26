function S = levint(Sn,k,n)

%kmax = length(Sn) - 2;
num = 0;
den = 0;
for j=0:k,
    coef = (-1)^j*nchoosek(k,j)*((n+j)/(n+k)).^(k-1);
    num = num + coef*Sn(n+j)/(Sn(n+j+1)-Sn(n+j));
    den = den + coef/(Sn(n+j+1)-Sn(n+j));
end
S = num/den;
