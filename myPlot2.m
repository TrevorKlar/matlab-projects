function myPlot2(t,y1, y2, plot_title, varname1, varname2)
%myPlot2 Plots 2 variables versus an independent variable
%   Expectations:
%   t = 1x(length) double
%   y1 = 1x(length) double
%   y2 = 1x(length) double
%   title = string          %Explains the plot in a word or two
%   varname1 = string        %What to call the plotted line in the legend
%   varname2 = string        %What to call the plotted line in the legend

    figure
    hold on % plot on the same figure
    plot(t,y1) 
    plot(t,y2) 
    hold off
    xlim([0 t(end)]) % stop a weird space from appearing on the right side of the plot
    title(strcat(plot_title))
    legend({strcat(varname1), strcat(varname2)})     
end