function output = thetaMethod(f, end_time, y_initial, time_step, theta)
%thetaMethod Solves a 2x2 ODE via the theta-method.
%   With some effort, this function could be rewritten to allow y1 to yk
%   using varargin, once i figure out how to pass syms into a function. 

syms t y1 y2 %Variables of y_prime

if theta ~= 1 %skip this for Forward Euler method
    Jf = matlabFunction(jacobian(f, [y1, y2])) %Symbolically differentiates 
    % to find the Jacobian, then converts from a symbolic expression to a
    % matlab function so we don't have to use any more symbolic methods.
end




% Initialize variables
start_time = 0;
time = start_time;
t_out = [start_time];    % these are the output data. We will graph y1_out 
y1_out = [y_initial(1)]; % and y2_out versus t_out.
y2_out = [y_initial(2)];

% Rename these constants for legibility 
h = time_step; 
tf = end_time; 

% Do stuff
for n = 1:(tf/h) %n matches the column index of output which was last filled.
    time = time + h; %increment time
    t_out = [t_out, time]; % append to time output
    g = [y1;y2] - [y1_out(n);y2_out(n)] - theta*h*subs(f, [t,y1,y2], [time-h,y1_out(n),y2_out(n)]) - (1-theta)*h*subs(f, t, time);
    % This would run waaaaaaay faster if I took advantage of the fact that
    % the Jacobian is always the same constant matrix for these problems, 
    % and precomputed it. That combined with letting g=matlabFunction(g) and
    % J=matlabFunction(J) and passing those to NewtonsMethod would huuugely
    % speed things up. Turns out, subs() is reeeeally slow, especially
    % doing it 1000 times. 
    new_y = NewtonsMethod2(g, [y1_out(n);y2_out(n)]); %compute next y
    y1_out = [y1_out, new_y(1)]; % append y to outputs
    y2_out = [y2_out, new_y(2)];
end

output = [t_out; y1_out; y2_out] ;
end

