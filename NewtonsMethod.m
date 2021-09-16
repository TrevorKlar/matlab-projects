function output = NewtonsMethod(F, J, root_guess)
%NewtonsMethod Finds a root of a function [F1(y1, y2); F2(y1, y2)] via Newton's Method
%   we're expecting F and J to be matlab functions of two variables here.
%   Here we require J to be computed in advance because we don't need a lot
%   of generality and precomputing is faster.

epsilon = 5*10^(-2); %How close to get before terminating

%Initializing variables
delta = [1;1]; % start with some arbitrary nonzero delta (this will change)
y = root_guess;
r = 0.5; % this initial value gives (r/(1-r))=1
n=0;

%do one Newton iteration no matter what
    %Plug in
    n=n+1;
    J_y = J(y(1), y(2));
    F_y = F(y(1), y(2));
    %Solve for delta in the equation (J_y)(delta)=-(F_y)
    delta0 = double(linsolve(J_y,-F_y));
    y = y + delta0; %save delta0 for use in future iterations
while (r/(1-r))*norm(delta)>epsilon 
    %Plug in
    n=n+1;
    J_y = J(y(1), y(2));
    F_y = F(y(1), y(2));
    %Solve for delta in the equation (J_y)(delta)=-(F_y)
    delta = double(linsolve(J_y,-F_y));
    y = y + delta;
    r = (norm(delta)/norm(delta0))^(1/n); %recalculate rate of convergence
end %loop until convergence test passes
output = y;
end

