function myPlot1(t,y, plot_title, varname)
%myPlot1 Plots 1 variable versus an independent variable
%   Expectations:
%   t = 1x(length) double
%   y = 1x(length) double
%   title = string          %Explains the plot in a word or two
%   varname = string        %What to call the plotted line in the legend

    figure
    plot(t,y) 
    xlim([0 t(end)]) % stop a weird space from appearing on the right side of the plot
    title(strcat(plot_title))
    legend({strcat(varname)})     
end