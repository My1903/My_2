function d_function = D(x,y)




dw  = 0.05;
dg  = 0.0025;
a = 1; 
b= 1;
R=1;
distance = x^2 + y^2;
dmask = distance <= R^2;
d_function= dmask .* a.*dw + (1-dmask) .* b.*dg;

% 
% if x.^2 + y.^2 >1
%     b =1;
% else 
%     b=0;
% end
% 
% d_function = a*dw*I + b*dg*I;

end







% a =piecewise(x<0.5 & y <0.5, 1, 0);
% b =piecewise(x>0.5 & y >0.5, 1, 0);
