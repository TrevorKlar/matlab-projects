% Theta method starts here

% Initialize variables
start_time = 0;
time = start_time;
t_out = [start_time];    % These are the output data. We will graph
y1_out = [y_initial(1)]; % y1_out and y2_out versus t_out.
y2_out = [y_initial(2)];
Jg = @(y1, y2) eye(2)  - (1-theta)*h*Jf;

% Compute values
for n = 1:(end_time/h) %n matches the column index of output which was last filled.
    time = time + h; %increment time
    t_out = [t_out, time]; % append to time output
    if theta ~= 1 % backward euler or trapezoidal method                
        g = @(y1,y2)( [y1;y2] - [y1_out(n);y2_out(n)] - theta*h*f(time-h,y1_out(n),y2_out(n)) - (1-theta)*h*f(time,y1,y2) );
        %Jg = @(y1, y2) eye(2)  - (1-theta)*h*Jf   
        % No need to run this line over and over, so we do it above 
        % the for block just once.

        % Now we compute next y values with Netwon's Method. 
        new_y = NewtonsMethod(g, Jg, [y1_out(n);y2_out(n)]); 

        y1_out = [y1_out, new_y(1)]; % append y to outputs
        y2_out = [y2_out, new_y(2)];
    else % forward Euler method
        new_y = [y1_out(n);y2_out(n)] + h*f(time-h,y1_out(n),y2_out(n));
        y1_out = [y1_out, new_y(1)]; % append y to outputs
        y2_out = [y2_out, new_y(2)];
    end
end

p = [t_out; y1_out; y2_out] ; % put all the outputs together for plotting
