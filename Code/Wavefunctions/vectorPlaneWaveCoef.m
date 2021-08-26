function [alm, blm] = vectorPlaneWaveCoef(L,Ex,Ey,Ez,thetak,phik)
% Multipole coefficients for plane wave at the origin
% for fully-normalized vector wave functions.
%
% L:                Maximum harmonic L
% Ex,Ey,Ez:         Electric field components
% thetak, phik:     k-vector direction (radians)
%
% alm, blm:         plane wave multipole coefficients, 
%                   harmonics 1:L, all m, linearly indexed
%
% Dependencies: sphericalY

tot = L^2 + 2*L;
alm = zeros(tot,1);
blm = zeros(tot,1);
khat = [sin(thetak)*cos(phik) sin(thetak)*sin(phik) cos(thetak)];
H = cross(khat,[Ex,Ey,Ez]);
Hx = H(1);
Hy = H(2);
Hz = H(3);
ylm = conj(sphericalY(L,thetak,phik));
e1 = 0.5*(Ex + 1i*Ey);
e2 = 0.5*(Ex - 1i*Ey);
h1 = 0.5*(Hx + 1i*Hy);
h2 = 0.5*(Hx - 1i*Hy);
for l=1:L,
    c1 = 4*pi*1i^(l+1)/sqrt(l*(l+1));
    for m=-l:l,
        elm = sqrt((l-m)*(l+m+1));
        flm = sqrt((l+m)*(l-m+1));
        ind = l^2 + l + m;
        atmp = Ez*m*ylm(ind);
        btmp = Hz*m*ylm(ind);
        if m < l
            atmp = atmp + e1*elm*ylm(ind+1);
            btmp = btmp + h1*elm*ylm(ind+1);
        end
        if m > -l
            atmp = atmp + e2*flm*ylm(ind-1);
            btmp = btmp + h2*flm*ylm(ind-1);
        end
        alm(ind) = c1*atmp;
        blm(ind) = 1i*c1*btmp;
    end
end