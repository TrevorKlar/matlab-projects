function output = myJacobian(f, vars)
%Jacobian find the Jacobian of a function

output = [];
for i = 1:length(vars)
   output = [output, diff(f, vars(i))];
end
