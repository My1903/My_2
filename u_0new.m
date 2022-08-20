%Initial condition
function u_0 = u_0new(X,Y)
 
%% Parameters setting
A = 2*10e2;
x0 =3;
y0 =3;
B = 0.6;
R0 = 1.8;
 
distance = (X - x0).^2 + (Y - y0).^2;
u0mask = distance <= R0^2;
u_0 = u0mask .* (A.* exp(-distance / B)) + (1-u0mask) .* 1;
%u_0new;
 
 
