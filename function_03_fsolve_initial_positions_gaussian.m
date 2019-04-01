function F = function_03_fsolve_initial_positions_gaussian( x )
%function_fsolve_initial_positions
%Determine the initial positions

global num_springs %total number of springs
global m_springs_per_cell %number of springs per cell
global L    %length of the domain
global gaussian_n %scaling of the distributon
global gaussian_sigma %standard deviation of the distribution
global gaussian_mu %mean of the distribution

F(1) = x(1) - 0;
for i=2:num_springs
    F(i) = (1/(x(i) - x(i-1)) + 1/(x(i) - x(i+1)))/(x(i-1)/2 - x(i+1)/2) + (m_springs_per_cell)*(gaussian_n*(x(i)-gaussian_mu)/(sqrt(2*pi)*gaussian_sigma^3))*(exp(-((x(i)-gaussian_mu)^2)/(2*gaussian_sigma^2))) ;
end
F(num_springs+1)= x(num_springs+1) - L;



end
