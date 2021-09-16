%% Project #1

syms t y1 y2 %Variables of y_prime

%% Problem 1
y_prime = [-y1;-100*(y2-sin(t))+cos(t)];
end_time = 1;
y_initial = [1;2];

for h= [0.05, 0.01]
    for theta = [0, 0.5, 1]
        p = thetaMethod(y_prime, end_time, y_initial, h, theta);

        figure
        hold on
        plot(p(1,:),p(2,:))
        plot(p(1,:),p(3,:))
        hold off
        xlim([0 end_time])
        title(strcat('Problem 1, h=', string(h), ', \theta=', string(theta)))
        legend({'y1', 'y2'})
    end
end
        
%% Problem 2
y_prime = [.25*y1-.01*y1*y2; -y2+.01*y1*y2];
end_time = 100;
y_initial = [10;10];

for h= [0.1, 0.001]
    for theta = [0, 0.5, 1]
        p = thetaMethod(y_prime, end_time, y_initial, h, theta);

        figure
        hold on
        plot(p(1,:),p(2,:))
        plot(p(1,:),p(3,:))
        hold off
        xlim([0 end_time])
        title(strcat('Problem 2, h=', string(h), ', \theta=', string(theta)))
        legend({'y1', 'y2'})
        
        figure
        hold on
        plot(p(2,:),p(3,:))
        hold off
        title(strcat('Problem 2, h=', string(h), ', \theta=', string(theta)))
        xlabel('y1')
        ylabel('y2')        
        
    end
end
%% Scratch
y_prime = [.25*y1-.01*y1*y2; -y2+.01*y1*y2];
end_time = 100;
y_initial = [10;10];

for h= [0.1, 0.001]
    for theta = [0, 0.5, 1]
        p = thetaMethod(y_prime, end_time, y_initial, h, theta);

        figure
        hold on
        plot(p(1,:),p(2,:))
        plot(p(1,:),p(3,:))
        hold off
        xlim([0 end_time])
        title(strcat('Problem 2, h=', string(h), ', \theta=', string(theta)))
        legend({'y1', 'y2'})
        
        figure
        hold on
        plot(p(2,:),p(3,:))
        hold off
        title(strcat('Problem 2, h=', string(h), ', \theta=', string(theta)))
        xlabel('y1')
        ylabel('y2')        
        
    end
end
