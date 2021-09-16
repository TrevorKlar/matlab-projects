%%%%%%% Project #2

% Options
plot_on = true;
ETOL = 1e-3;

%% Problem 1 - Easy Problem

% The given data
f = @(t,y)[-y(1);-10*(y(2)-t^2)+2*t];
end_time = 1;
y_initial = [1;2];

[t,y,h_list,f_list] = PECE_MethodOrder2(f, y_initial, end_time, ETOL);

offStepPoint(1, t, y, h_list, f_list)

if plot_on    
    myPlot2(t, y(1,:), y(2,:), 'Problem 1 PECE', 'y1', 'y2');
    exportgraphics(gcf, 'prob1_y1_and_y2_vs_t.png')
    myPlot1(t, h_list, 'Problem 1 PECE', 'h');     
    exportgraphics(gcf, 'prob1_h.png')
end



%% Problem 2 - Predator-Prey Problem

% The given data
f = @(t,y)[.25*y(1)-.01*y(1)*y(2); -y(2)+.01*y(1)*y(2)];
%Jf = @(y1,y2)[.25 - .01*y2, -.01*y1; .01*y2, .01*y1 - 1];
end_time = 100;
y_initial = [10;10];

[t,y,h_list,f_list] = PECE_MethodOrder2(f, y_initial, end_time, ETOL);

offStepPoint(1, t, y, h_list, f_list)

if plot_on    
    myPlot2(t, y(1,:), y(2,:), 'Problem 2 PECE', 'y1', 'y2')
    xlabel('t')
    exportgraphics(gcf, 'prob2_y1_and_y2_vs_t.png')
    figure
        plot(y(1,:),y(2,:)) %y2 vs y1
        title('Problem 2 PECE')
        xlabel('y1')
        ylabel('y2')   
        exportgraphics(gcf, 'prob2_y1_vs_y2.png')
    myPlot1(t, h_list, 'Problem 2 PECE stepsizes', 'h');     
    exportgraphics(gcf, 'prob2_h.png')
end

%% Problem 3 - Van der Pol's Equation

% The given data
eta = 2
f = @(t,y)[y(2);eta*((1-y(1)^2)*y(2)-y(1))];
end_time = 11;
y_initial = [2;0];

[t,y,h_list,f_list] = PECE_MethodOrder2(f, y_initial, end_time, ETOL);

offStepPoint(1, t, y, h_list, f_list)

if plot_on    
    myPlot2(t, y(1,:), y(2,:), 'Problem 3 PECE', 'y1', 'y2')
    xlabel('t')
    exportgraphics(gcf, 'prob3_y1_and_y2_vs_t.png')
    figure
        plot(y(1,:),y(2,:)) %y2 vs y1
        title('Problem 3 PECE')
        xlabel('y1')
        ylabel('y2')   
        exportgraphics(gcf, 'prob3_y1_vs_y2.png')
    myPlot1(t, h_list, 'Problem 3 PECE stepsizes', 'h');     
    exportgraphics(gcf, 'prob3_h.png')
end

%% Problem 4 - Advection PDE via MOL

% Knobs I can turn
partition_count=100;

% Using MOL on the given data to construct an ODE
spatial_domain = [0,1];
delta_x=(spatial_domain(2)-spatial_domain(1))/partition_count;
%Note that u_t + u_x = 0 for this problem, so u_t = -u_x, which we will
%approximate. 

%Now we construct the coefficient matrix using the first order backward
%difference approximation to the derivative u_x.
A = zeros(partition_count);
    % start at i=2 since u_0=1 always, so u_0'=0
    for i = 2:partition_count %u_i' = u_{i-1}-u_i, over delta_x (done later)
        A(i, i-1) = 1;
        A(i, i) = -1;
    end
f = @(t, y)[(1/delta_x)*A*y]; % here we divide by delta_x
end_time = 1;
x = linspace(spatial_domain(1), spatial_domain(2), partition_count)'; 
y_initial = exp(-10*x);

times = [0,0.25,0.5,0.6,0.8,1]; %times at which we want to find u(x,t) for 
                                %each t in the list
                                
[t,y,h_list,f_list] = PECE_MethodOrder2(f, y_initial, end_time, ETOL);

y_offStepPoints = zeros(partition_count, length(times));
for i=1:length(times)
    y_offStepPoints(:,i) = offStepPoint(times(i), t, y, h_list, f_list);
end

%plot u vs t with lines ranging across x
figure
hold on % plot on the same figure
cc=parula(partition_count);
for i=1:partition_count
    plot(t, y(i,:), 'color', cc(i,:))
end
c = colorbar('Direction','reverse');
c.Label.String = "color denotes x coordinate";
hold off
xlim([0 t(end)]) % stop a weird space from appearing on the right side of the plot
ylim([0 1])
title(strcat('Problem 4 MOL'))
xlabel('t, time')
ylabel('u, temperature')
exportgraphics(gcf, 'prob4_u_vs_t.png')
%legend({strcat(varname1), strcat(varname2)})     

%plot u vs x with lines ranging across t
figure
hold on % plot on the same figure
cc=copper(length(times));
for i=1:length(times)
    plot(x, y_offStepPoints(:,i), 'color', cc(i,:))
end
hold off
xlim([0 times(end)]) % stop a weird space from appearing on the right side of the plot
ylim([0 1])
title(strcat('Problem 4 MOL, Offstep Points'))
xlabel('x')
ylabel('u')
legend({strcat("t=", string(times(1))), strcat("t=", string(times(2))), ...
    strcat("t=", string(times(3))), strcat("t=", string(times(4))), ...
    strcat("t=", string(times(5))), strcat("t=", string(times(6))), })
exportgraphics(gcf, 'prob4_u_vs_x.png')

%plot h
myPlot1(t, h_list, 'Problem 4 MOL stepsizes', 'h');   
xlabel('t')
ylabel('h')
exportgraphics(gcf, 'prob4_h.png')