

%%hypothetical values
%u_t = -u_x;


%boundary_value = 1                     %u(t,0)=1
%initial_value_function = exp(-10x)     %u(0,x)=exp(-10x)

%%intialize variables
spatial_domain = [0,1];
partition_count=10;
delta_x=(spatial_domain(2)-spatial_domain(1))/partition_count;
A = zeros(partition_count);
    for i = 2:partition_count
        A(i, i-1) = 1;
        A(i, i) = -1;
    end
f = @(t, y)[(1/delta_x)*A*y];
end_time = 1;
x = linspace(spatial_domain(1), spatial_domain(2), partition_count); 
y_initial = exp(-10*x);


