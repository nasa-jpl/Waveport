


%%  

syms r_x r_y y_1 x_2 y_2 n_1 n_2



leftnum = r_x^2*y_1^2;
leftden = r_x^2+(y_1-r_y)^2;
rightnum = n_2^2*(r_x*y_2-r_y*x_2)^2;
rightden = n_1^2*((x_2-r_x)^2+(y_2-r_y)^2);

M = leftnum*rightden - rightnum*leftden;

collect(simplify(M),'r_x')


r = [r_x r_y 0];
r1 = [0 y_1 0];
r2 = [x_2 y_2 0];
v1 = r1 - r;
v2 = r2 - r;
sint1 = cross(r,v1)/norm(r)/norm(v1);
sint2 = cross(-r,v2)/norm(r)/norm(v2);
sint1 = sint1(3)^2;
sint2 = sint2(3)^2;
 
leftnum = (r_x*r_y - r_x*(r_y - y_1))^2;
leftden = (abs(r_x)^2 + abs(r_y)^2)*(abs(r_y - y_1)^2 + abs(r_x)^2);
rightnum = (r_y*(r_x - x_2) - r_x*(r_y - y_2))^2;
rightden = (abs(r_x - x_2)^2 + abs(r_y - y_2)^2)*(abs(r_x)^2 + abs(r_y)^2);

leftnum = (r_x*r_y - r_x*(r_y - y_1))^2;
leftden = (r_y - y_1)^2 + r_x^2;
rightnum = (r_y*(r_x - x_2) - r_x*(r_y - y_2))^2;
rightden = (r_x - x_2)^2 + (r_y - y_2)^2;

r_y = sqrt(r^2 - r_x^2);
leftnum = r_x^2*y_1^2;
leftden = (r_y - y_1)^2 + r_x^2;
rightnum = (r_x*y_2 - r_y*x_2)^2;
rightden = (r_x - x_2)^2 + (r_y - y_2)^2;


M = n_1^2*leftnum*rightden - n_2^2*rightnum*leftden;

T = collect(simplify(expand(M)),'r_x')

clear r
syms r
r_y = sqrt(r^2 - r_x^2)


%%

syms c_1 c_2 c_3 c_4 c_5 c_6 c_7 c_8 c_9 x y r

leftnum = c_1*x^2;
leftden = c_2 + c_3*y;
rightnum = c_4*x^2 + c_5*x*y + c_6;
rightden = c_7 + c_8*x + c_9*y;

N = leftnum*rightden - rightnum*leftden
collect(expand(N),'y')

N = (-c_3*c_5*x)*y^2 + (c_1*c_9*x^2 - c_2*c_5*x - c_3*c_4*x^2 - c_3*c_6)*y + c_1*c_7*x^2 - c_2*c_4*x^2 - c_2*c_6 + c_1*c_8*x^3

% sub y^2
M2 = (-c_3*c_5*x)*(r^2-x^2) + (c_1*c_9*x^2 - c_2*c_5*x - c_3*c_4*x^2 - c_3*c_6)*y + c_1*c_7*x^2 - c_2*c_4*x^2 - c_2*c_6 + c_1*c_8*x^3
% expand y^2 term
M2 = c_1*c_7*x^2 - y*(c_3*c_6 + c_2*c_5*x + c_3*c_4*x^2 - c_1*c_9*x^2) - c_2*c_4*x^2 - c_2*c_6 + c_1*c_8*x^3 - c_3*c_5*x*(r^2 - x^2)
collect(expand(M2),'y')
(c_1*c_9*x^2 - c_2*c_5*x - c_3*c_4*x^2 - c_3*c_6)*y + c_1*c_7*x^2 - c_2*c_4*x^2 - c_2*c_6 + c_3*c_5*x^3 + c_1*c_8*x^3 - c_3*c_5*r^2*x


Ma = Mb*y;

Ma = c_1*c_7*x^2 - c_2*c_4*x^2 - c_2*c_6 + c_3*c_5*x^3 + c_1*c_8*x^3 - c_3*c_5*r^2*x
Mb = -(c_1*c_9*x^2 - c_2*c_5*x - c_3*c_4*x^2 - c_3*c_6)

collect(Ma,'x')
(c_3*c_5 + c_1*c_8)*x^3 + (c_1*c_7 - c_2*c_4)*x^2 - c_3*c_5*r^2*x  - c_2*c_6

collect(Mb)
(c_3*c_4 - c_1*c_9)*x^2 + c_2*c_5*x + c_3*c_6

b_1 = c_3*c_5 + c_1*c_8
b_2 = c_1*c_7 - c_2*c_4
b_3 = -c_3*c_5*r^2
b_4 = -c_2*c_6
b_5 = c_3*c_4 - c_1*c_9
b_6 = c_2*c_5
b_7 = c_3*c_6


syms b_1 b_2 b_3 b_4 b_5 b_6 b_7
left = b_1*x^3 + b_2*x^2  + b_3*x + b_4;
right = (b_5*x^2 + b_6*x + b_7)*y;

O = left^2 - right^2
O = (b_1*x^3 + b_2*x^2 + b_3*x + b_4)^2 - y^2*(b_5*x^2 + b_6*x + b_7)^2
P = (b_1*x^3 + b_2*x^2 + b_3*x + b_4)^2 - (r^2-x^2)*(b_5*x^2 + b_6*x + b_7)^2

collect(expand(P),'x')


(b_1^2 + b_5^2)*x^6 + (2*b_1*b_2 + 2*b_5*b_6)*x^5 + (b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3)*x^4 + (- 2*b_5*b_6*r^2 + 2*b_1*b_4 + 2*b_2*b_3 + 2*b_6*b_7)*x^3 + (b_3^2 - b_6^2*r^2 + b_7^2 - 2*b_5*b_7*r^2 + 2*b_2*b_4)*x^2 + (- 2*b_6*b_7*r^2 + 2*b_3*b_4)*x + b_4^2 - b_7^2*r^2

b_1^2 + b_5^2 
2*b_1*b_2 + 2*b_5*b_6
b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3
-2*b_5*b_6*r^2 + 2*b_1*b_4 + 2*b_2*b_3 + 2*b_6*b_7 
b_3^2 - b_6^2*r^2 + b_7^2 - 2*b_5*b_7*r^2 + 2*b_2*b_4
-2*b_6*b_7*r^2 + 2*b_3*b_4
b_4^2 - b_7^2*r^2


a_6*x^6 + a_5*x^5 + a_4*x^4 + a_3*x^3 + a_2*x^2 + a_1*x + a_0


a_6 = b_1^2 + b_5^2 
a_5 = 2*b_1*b_2 + 2*b_5*b_6
a_4 = b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3
a_3 = -2*b_5*b_6*r^2 + 2*b_1*b_4 + 2*b_2*b_3 + 2*b_6*b_7 
a_2 = b_3^2 - b_6^2*r^2 + b_7^2 - 2*b_5*b_7*r^2 + 2*b_2*b_4
a_1 = -2*b_6*b_7*r^2 + 2*b_3*b_4
a_0 = b_4^2 - b_7^2*r^2





%%

r = 1;
h = 1;
n1 = 1;
n2 = 2;

y1 = r+h;
x2 = 0.4;
y2 = 0.4;

figure(1),clf,hold all
scatter(0,y1,'filled')
scatter(x2,y2,'filled')
plotcircle(0,0,r)
hold off
axis equal

ry = r;
for n=1:10,

a = n1.^2.*y1.^2 - n2.^2.*y2.^2;
b = -2*n1.^2.*x2.*y1.^2 + 2*n2.^2.*ry.*x2.*y2;
c = n1.^2.*y1.^2.*((ry-y2).^2 + x2.^2) - n2.^2.*y2.^2.*(ry-y1).^2 - n2.^2.*ry.^2.*x2.^2;
d = 2*n2.^2.*ry.*x2.*y2.*(ry-y1).^2;
e = n2.^2.*ry.^2.*x2.^2.*(ry-y1).^2;

[w1 w2 w3 w4] = quarticroots(a,b,c,d,e);
r = [w1(:) w2(:) w3(:) w4(:)];
u = zeros(N,1);
for n=1:4,
    ind = and(and(imag(r(:,n))==0,r(:,n)>0),r(:,n)<L(:));
    u(ind) = r(ind,n);
end
u = reshape(u,sz);
end


%%



%%

clear

r = 1;
h = 1;
n1 = 1;
n2 = 2;

x1 = -2;
y1 = 1;
x2 = -0.5;
y2 = 0.4;

figure(1),clf,hold all
scatter(0,y1,'filled')
scatter(x2,y2,'filled')
plotcircle(0,0,r)
hold off
axis equal


[rx ry theta1 theta2 v1mag v2mag] = refractionCircle(x1,y1,x2,y2,n1,n2,r);

figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
scatter([rx],[ry],'filled')
scatter(0,0,'filled')
plotcircle(0,0,r)
%plot([0 0],[0 y1])
plot([x2 rx],[y2 ry])
plot([x1 rx],[y1 ry])
%plot([0 rx],[0 ry])
hold off
axis equal


%% version 2

clear
syms r_x r_y x_1 y_1 x_2 y_2 n_1 n_2 r G

rr = [r_x r_y];
r1 = [x_1 y_1];
r2 = [x_2 y_2];
v1 = r1 - rr;
v2 = r2 - rr;

v1 = sqrt(sum(v1.^2));
v2 = sqrt(sum(v2.^2));

r_y = sqrt(r^2 - r_x^2);
L = n_1*v1 + n_2*v2;
L = eval(L)
diff(L,'r_x')


left = (n_1*(2*r_x - 2*x_1 + (2*r_x*(y_1 - (r^2 - r_x^2)^(1/2)))/(r^2 - r_x^2)^(1/2)))/(2*((r_x - x_1)^2 + (y_1 - (r^2 - r_x^2)^(1/2))^2)^(1/2))
right = -(n_2*(2*r_x - 2*x_2 + (2*r_x*(y_2 - (r^2 - r_x^2)^(1/2)))/(r^2 - r_x^2)^(1/2)))/(2*((r_x - x_2)^2 + (y_2 - (r^2 - r_x^2)^(1/2))^2)^(1/2))


M1 = simplify(left^2);
M2 = simplify(right^2);

M1 = (n_1^2*(2*r_x - 2*x_1 + (2*r_x*(y_1 - (r^2 - r_x^2)^(1/2)))/(r^2 - r_x^2)^(1/2))^2)/(4*((r_x - x_1)^2 + (y_1 - (r^2 - r_x^2)^(1/2))^2));
M2 = (n_2^2*(2*r_x - 2*x_2 + (2*r_x*(y_2 - (r^2 - r_x^2)^(1/2)))/(r^2 - r_x^2)^(1/2))^2)/(4*((r_x - x_2)^2 + (y_2 - (r^2 - r_x^2)^(1/2))^2));

N1 = simplify((r^2-r_x^2)^2*M1)/(r^2-r_x^2);
N2 = simplify((r^2-r_x^2)^2*M2)/(r^2-r_x^2);

O1 = simplify(expand(N1))
O2 = simplify(expand(N2))

collect(O1,'r_x')
collect(O2,'r_x')

O1 = ((n_1^2*x_1^2 - n_1^2*y_1^2)*r_x^2 + 2*n_1^2*x_1*y_1*(r^2 - r_x^2)^(1/2)*r_x - n_1^2*r^2*x_1^2)/(2*x_1*r_x + 2*y_1*(r^2 - r_x^2)^(1/2) - r^2 - x_1^2 - y_1^2)
O2 = ((n_2^2*x_2^2 - n_2^2*y_2^2)*r_x^2 + 2*n_2^2*x_2*y_2*(r^2 - r_x^2)^(1/2)*r_x - n_2^2*r^2*x_2^2)/(2*x_2*r_x + 2*y_2*(r^2 - r_x^2)^(1/2) - r^2 - x_2^2 - y_2^2)

O1 = (c_1*r_x^2 + c_2*(r^2 - r_x^2)^(1/2)*r_x + c_3)/(c_4*r_x + c_5*(r^2 - r_x^2)^(1/2) + c_6)
O2 = (c_7*r_x^2 + c_8*(r^2 - r_x^2)^(1/2)*r_x + c_9)/(c_10*r_x + c_11*(r^2 - r_x^2)^(1/2) + c_12)

c_1 = (n_1^2*x_1^2 - n_1^2*y_1^2);
c_2 = 2*n_1^2*x_1*y_1
c_3 = -n_1^2*r^2*x_1^2
c_4 = 2*x_1
c_5 = 2*y_1
c_6 = -r^2 - x_1^2 - y_1^2
c_7..c_12 same as c_1..c_6 except 1->2

syms c_1 c_2 c_3 c_4 c_5 c_6 c_7 c_8 c_9 c_10 c_11 c_12  g

O1 = (c_1*r_x^2 + c_2*(r^2 - r_x^2)^(1/2)*r_x + c_3)/(c_4*r_x + c_5*(r^2 - r_x^2)^(1/2) + c_6)
O2 = (c_7*r_x^2 + c_8*(r^2 - r_x^2)^(1/2)*r_x + c_9)/(c_10*r_x + c_11*(r^2 - r_x^2)^(1/2) + c_12)


leftnum = (c_1*r_x^2 + c_2*(r^2 - r_x^2)^(1/2)*r_x + c_3);
leftden = (c_4*r_x + c_5*(r^2 - r_x^2)^(1/2) + c_6);
rightnum = (c_7*r_x^2 + c_8*(r^2 - r_x^2)^(1/2)*r_x + c_9);
rightden = (c_10*r_x + c_11*(r^2 - r_x^2)^(1/2) + c_12);

leftnum = (c_1*r_x^2 + c_2*g*r_x + c_3);
leftden = (c_4*r_x + c_5*g + c_6);
rightnum = (c_7*r_x^2 + c_8*g*r_x + c_9);
rightden = (c_10*r_x + c_11*g + c_12);

P1 = leftnum*rightden
P2 = rightnum*leftden

collect(simplify(expand(P1)),'g')
collect(simplify(expand(P2)),'g')

Q1 = c_2*c_11*r_x*(r^2 - r_x^2) + (c_11*(c_1*r_x^2 + c_3) + c_2*r_x*(c_12 + c_10*r_x))*g + (c_1*r_x^2 + c_3)*(c_12 + c_10*r_x)
Q2 = c_5*c_8*r_x*(r^2 - r_x^2) + (c_5*(c_7*r_x^2 + c_9) + c_8*r_x*(c_6 + c_4*r_x))*g + (c_7*r_x^2 + c_9)*(c_6 + c_4*r_x)

collect(simplify(expand(Q1)),'g')
collect(simplify(expand(Q2)),'g')

(c_3*c_11 + c_2*c_12*r_x + c_1*c_11*r_x^2 + c_2*c_10*r_x^2)*g + c_3*c_12 + c_3*c_10*r_x + c_1*c_10*r_x^3 + c_1*c_12*r_x^2 - c_2*c_11*r_x^3 + c_2*c_11*r^2*r_x
(c_5*c_9 + c_6*c_8*r_x + c_4*c_8*r_x^2 + c_5*c_7*r_x^2)*g + c_6*c_9 + c_4*c_9*r_x + c_4*c_7*r_x^3 + c_6*c_7*r_x^2 - c_5*c_8*r_x^3 + c_5*c_8*r^2*r_x


% collect r_x in g
((c_1*c_11 + c_2*c_10)*r_x^2 + c_2*c_12*r_x + c_3*c_11)*g  + (c_1*c_10 - c_2*c_11)*r_x^3 + c_1*c_12*r_x^2  + (c_3*c_10 + c_2*c_11*r^2)*r_x + c_3*c_12 
((c_4*c_8 + c_5*c_7)*r_x^2 + c_6*c_8*r_x + c_5*c_9)*g + (c_4*c_7 - c_5*c_8)*r_x^3 + c_6*c_7*r_x^2 + (c_4*c_9 + c_5*c_8*r^2)*r_x + c_6*c_9  


% put g on the right side
left3 = (c_1*c_10 - c_2*c_11)*r_x^3 + c_1*c_12*r_x^2  + (c_3*c_10 + c_2*c_11*r^2)*r_x + c_3*c_12 - ((c_4*c_7 - c_5*c_8)*r_x^3 + c_6*c_7*r_x^2 + (c_4*c_9 + c_5*c_8*r^2)*r_x + c_6*c_9)
right3 = ((c_4*c_8 + c_5*c_7)*r_x^2 + c_6*c_8*r_x + c_5*c_9)*g - ((c_1*c_11 + c_2*c_10)*r_x^2 + c_2*c_12*r_x + c_3*c_11)*g 

collect(left3,'r_x')
collect(right3,'g')

(c_1*c_10 - c_4*c_7 - c_2*c_11 + c_5*c_8)*r_x^3 + (c_1*c_12 - c_6*c_7)*r_x^2 + (c_3*c_10 - c_4*c_9 + c_2*c_11*r^2 - c_5*c_8*r^2)*r_x + c_3*c_12 - c_6*c_9

(c_5*c_9 - c_3*c_11 - r_x^2*(c_1*c_11 + c_2*c_10) + r_x^2*(c_4*c_8 + c_5*c_7) - c_2*c_12*r_x + c_6*c_8*r_x)*g
((c_4*c_8 + c_5*c_7)-(c_1*c_11 + c_2*c_10))*r_x^2  + (c_6*c_8 - c_2*c_12)*r_x + c_5*c_9 - c_3*c_11  )*g

% b coefficients
b_1*r_x^3 + b_2*r_x^2 + b_3*r_x + b_4
(b_5*r_x^2  + b_6*r_x + b_7)*g


b_1 = c_1*c_10 - c_4*c_7 - c_2*c_11 + c_5*c_8
b_2 = c_1*c_12 - c_6*c_7
b_3 = c_3*c_10 - c_4*c_9 + (c_2*c_11 - c_5*c_8)*r^2
b_4 = c_3*c_12 - c_6*c_9
b_5 = c_4*c_8 + c_5*c_7 - c_1*c_11 - c_2*c_10
b_6 = c_6*c_8 - c_2*c_12
b_7 = c_5*c_9 - c_3*c_11  


clear 
syms b_1 b_2 b_3 b_4 b_5 b_6 b_7 r r_x
% b coefficients
left4 = b_1*r_x^3 + b_2*r_x^2 + b_3*r_x + b_4;
right4 = (b_5*r_x^2  + b_6*r_x + b_7)*sqrt(r^2 - r_x^2);

% square 

collect(left4^2 - right4^2,'r_x')

(b_1^2 + b_5^2)*r_x^6 + (2*b_1*b_2 + 2*b_5*b_6)*r_x^5 + (b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3)*r_x^4 + (- 2*b_5*b_6*r^2 + 2*b_1*b_4 + 2*b_2*b_3 + 2*b_6*b_7)*r_x^3 + (2*b_2*b_4 - r^2*(b_6^2 + 2*b_5*b_7) + b_3^2 + b_7^2)*r_x^2 + (- 2*b_6*b_7*r^2 + 2*b_3*b_4)*r_x + b_4^2 - b_7^2*r^2

a_6 = b_1^2 + b_5^2
a_5 = 2*(b_1*b_2 + b_5*b_6)
a_4 = b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3
a_3 = 2*(-b_5*b_6*r^2 + b_1*b_4 + b_2*b_3 + b_6*b_7)
a_2 = 2*b_2*b_4 - r^2*(b_6^2 + 2*b_5*b_7) + b_3^2 + b_7^2
a_1 = 2*(-b_6*b_7*r^2 + b_3*b_4)
a_0 = b_4^2 - b_7^2*r^2



%%  Iterative











