function [y y_real] = data_generating(x, type)

a = 3.14;
aa = 0.5;
randprob = 0.7;

if type == "sin"
   p = rand();
   if p > randprob
   y = sin(x) * p;
   else
   y = sin(x);
   end
   y_real = sin(x);
elseif type == "cos"
   p = rand();
   if p > randprob
       y = cos(x) * p;
   else
       y = cos(x);
   end
   y_real = cos(x);
elseif type == "linear"
   p = rand();
   if p > randprob
       y = a * x * p;
   else
       y = a * x;
   end
   y_real = a * x;
elseif type == "nonlinear"
   p = rand();
   if p > randprob
       y = (x^2/aa + x^2/(1-aa)) * p;
   else
       y = x^2/aa + x^2/(1-aa);
   end
   y_real = x^2/aa + x^2/(1-aa);
end

end