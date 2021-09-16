%%%%%%% Project #1

% Options
plot_on = true;

%% Problem 1

% The given data
f = @(t,y1,y2)[-y1;-100*(y2-sin(t))+cos(t)];
Jf = [-1,0;0,-100]; % Jacobian of f. No need to have the software compute 
                    % it since it's constant. 
end_time = 1;
y_initial = [1;2];

for h= [0.05, 0.01]
    for theta = [0, 0.5, 1]
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
        
        if plot_on
            figure
            hold on % plot on the same figure
            plot(p(1,:),p(2,:)) % y1 vs t
            plot(p(1,:),p(3,:)) % y2 vs t
            hold off
            xlim([0 end_time]) % stop a weird space from appearing on the right side of the plot
            title(strcat('Problem 1, h=', string(h), ', \theta=', string(theta)))
            legend({'y1', 'y2'})
        end
    end
end


           
        
%% Problem 2

% The given data
f = @(y1,y2)[.25*y1-.01*y1*y2; -y2+.01*y1*y2];
Jf = @(y1,y2)[.25 - .01*y2, -.01*y1; .01*y2, .01*y1 - 1];
end_time = 100;
y_initial = [10;10];

for h= [0.1, 0.001]
    for theta = [0, 0.5, 1]
        % Theta method starts here
        
        % Initialize variables
        start_time = 0;
        time = start_time;
        t_out = [start_time];    % These are the output data. We will graph y1_out 
        y1_out = [y_initial(1)]; % and y2_out versus t_out.
        y2_out = [y_initial(2)];
        Jg = @(y1, y2) eye(2)  - (1-theta)*h*Jf(y1,y2);
              
        % Compute values
        for n = 1:(end_time/h) %n matches the column index of output which was last filled.
            time = time + h; %increment time
            t_out = [t_out, time]; % append to time output
            if theta ~= 1 % backward euler or trapezoidal method   
                g = @(y1,y2)( [y1;y2] - [y1_out(n);y2_out(n)] - theta*h*f(y1_out(n),y2_out(n)) - (1-theta)*h*f(y1,y2) );
                %Jg = @(y1, y2) eye(2)  - (1-theta)*h*Jf(y1,y2)
                % No need to run this line over and over, so we do it above 
                % the for block just once.

                % Now we compute next y values with Netwon's Method. 
                new_y = NewtonsMethod(g, Jg, [y1_out(n);y2_out(n)]);
                y1_out = [y1_out, new_y(1)]; % append y to outputs
                y2_out = [y2_out, new_y(2)];
            else % forward Euler method
                new_y = [y1_out(n);y2_out(n)] + h*f(y1_out(n),y2_out(n));
                y1_out = [y1_out, new_y(1)]; % append y to outputs
                y2_out = [y2_out, new_y(2)];
            end %if
        end %for n

        p = [t_out; y1_out; y2_out] ; % put all the outputs together for plotting

        if plot_on    
            figure
            hold on % plot on the same figure
            plot(p(1,:),p(2,:)) % y1 vs t
            plot(p(1,:),p(3,:)) % y2 vs t
            hold off
            xlim([0 end_time]) % stop a weird space from appearing on the right side of the plot
            title(strcat('Problem 2, h=', string(h), ', \theta=', string(theta)))
            legend({'y1', 'y2'})

            figure
            %hold on
            plot(p(2,:),p(3,:)) %y2 vs y1
            %hold off
            title(strcat('Problem 2, h=', string(h), ', \theta=', string(theta)))
            xlabel('y1')
            ylabel('y2')   
        end %if
    end %for theta
end %for h

function output = NewtonsMethod(F, J, root_guess)
%NewtonsMethod Finds a root of a function [F1(y1, y2); F2(y1, y2)] via Newton's Method
%   we're expecting F and J to be matlab functions of two variables here.
%   Here we require J to be computed in advance because we don't need a lot
%   of generality and precomputing is faster.

epsilon = 5*10^(-2); %How close to get before terminating

%Initializing variables
delta = [1;1]; % start with some arbitrary nonzero delta (this will change)
y = root_guess;

while norm(delta)>epsilon
    %Plug in
    J_y = J(y(1), y(2));
    F_y = F(y(1), y(2));
    %Solve for delta in the equation (J_y)(delta)=-(F_y)
    delta = double(linsolve(J_y,-F_y));
    y = y + delta;
end %loop until delta <= epsilon
output = y;
end

