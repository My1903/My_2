% function f= f(u)
% f= log(abs(u));
% end
function f= f(u)
  EPSILON = 1e-6;
%   u = abs(u);
  if u < EPSILON
    f = 0;
  else 
   
     pro=1;
    f = pro*u*log(u);
  
  end
end
