function [t_out, y_out, h_out, f_out] = PECE_MethodOrder2(f, y_initial, end_time, error_tolerance)
%thetaMethod Solves an ODE via the 2nd order PECE method, with variable 
%   time steps.
%
%   expectations:
%   f   @(t, y)function 
%       y   Nx1 double
%       t   double

% Initialize variables
start_time = 0;
t_out = [start_time];    % These are the output data. We will graph
y_out = [y_initial]; %     y1_out and y2_out versus t_out.
f_out = [f(start_time, y_initial)];
h_out = [0];
AB_C3 = 5/12; % Coefficients for Adams-Bashforth and Adams-Moulton methods
AM_C3 = -1/12;%     respectively, taken from tables 5.1 and 5.2 in the text.
C = AM_C3/(AM_C3-AB_C3);

% Compute 2nd values
h = .0001; % We need another "intial" value to get PECE going
%FEM with tiny h
    new_time = t_out(end) + h; 
    new_y = y_out(:,end) + h*f_out(:,end);
    new_f = f(new_time, new_y); 
    t_out = [t_out, new_time]; % append to outputs
    y_out = [y_out, new_y]; 
    f_out = [f_out, new_f];
    h_out = [h_out, h];
%end FEM
    
% Compute remaining values
while t_out(end) < end_time %n matches the column index of output which was last filled.
    new_time = t_out(end) + h; %increment time
    new_y_p = y_out(:,end) + (h/2)*(3*f_out(:,end)-f_out(:,end-1)); %Predict via AB
    new_f = f(new_time, new_y_p); %Evaluate f
    new_y = y_out(:,end) + (h/2)*f_out(:,end) + (h/2)*new_f; %Correct via AM
    new_f = f(new_time, new_y); %Evaluate f
    error_estimate = C*norm(new_y - new_y_p);
    if false %error_estimate > error_tolerance 
        error_estimate
    end
    if error_estimate < error_tolerance 
        t_out = [t_out, new_time]; %append to outputs
        y_out = [y_out, new_y]; 
        f_out = [f_out, new_f];
        h_out = [h_out, h];
    end 
    h = computeTimestep(h, error_estimate, error_tolerance); 
end
h;
end

%% Function definitions
function output = computeTimestep(h, error_estimate, error_tolerance)
%computeTimestep Computes the size of the next timestep for PECE.
%   This function returns 
%        h*(tolerance_margin*error_tolerance/error_estimate)^(1/3). 
%
%   We want to choose h so that our error estimate is within our error
%   tolerance, but close to it (say, 90%). Thus we set 
%                       tolerance_margin = 0.9
%   We take the cube root because this is a 2nd order method which means
%   the error is of order h^3, thus our r-value is going to get cubed. 

tolerance_margin = 0.9; %default 0.9

r = (tolerance_margin*error_tolerance/error_estimate)^(1/3);

output=r*h;
end



