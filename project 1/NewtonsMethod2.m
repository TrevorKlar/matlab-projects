function output = NewtonsMethod2(F, root_guess)
%NewtonsMethod Finds a root of a function [F1(y1, y2); F2(y1, y2)] via Newton's Method

syms y1 y2 %Variables of f
epsilon = 5*10^(-2); %How close to get before terminating

%Initializing variables
delta = [1;1]; % start with some arbitrary nonzero delta (this will change)
y = root_guess;

%Take Jacobian 
J = myJacobian(F, [y1, y2]);

while norm(delta)>epsilon
    %Plug in
    J_y = double(subs(J, [y1, y2], [y(1), y(2)]));
    F_y = double(subs(F, [y1, y2], [y(1), y(2)]));
    %Solve for delta in the equation (J_y)(delta)=-(F_y)
    delta = double(linsolve(J_y,-F_y));
    y = y + delta;
end %loop until delta <= epsilon
output = y;
end

