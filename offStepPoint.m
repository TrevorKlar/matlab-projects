function output = offStepPoint(t_star,t,y,h,f)
%offStepPoint Interpolates to find y(t) for any valid t, given the data
%output from PECE_MethodOrder2(). 
%   Expectations: 
%   t = double, with t_list(1)<t<t_list(end). This is the point at which we
%       are evaluating.

n=find(t<=t_star); %find indeces of t_list for elements <= t
n=n(end); %only keep the largest such index
n=n+1; %in the theory we're used to thinking of y(n) as the value we're 
        %about to compute, and y(n-k) as past values, so now the last
        %computed value is y(n-1). 
if t_star==t(n-1) %check to see if we're done already
    output = y(:,n-1);
    return
end
h_star = t_star-t(n-1);
% This formula is obtained by approximating f at t(n-1) and t(n-2) by the
% interpolating polynomial, and integrating both sides from t(n-1) to
% t_star. Below we simply plug in to the formula obtained on paper. 
y_star = y(:,n-1)+h_star*f(:,n-2)+((f(:,n-2)-f(:,n-1))/(2*h(n-1)))*h_star^2;
output=y_star;
end


