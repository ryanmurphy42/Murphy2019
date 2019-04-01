function function_05_pde_spring_q_k_a( filepath_save ,time_vector_record )
%function_05_pde_spring_q_k_a - simulate the continuum model.

%% Initial conditions

%Density condition - normal distribution conditions
global gaussian_n
global gaussian_sigma
global gaussian_mu

gaussian_n=11.056818502502322;
gaussian_sigma=sqrt(9);
gaussian_mu=5;

%Spatial step
dx=0.05;

%Timestep initial
dt_ori=0.001;
dt_courant=0.0001;

%Timestep to ensure record near timepoints
dt_near_time_record = 0.0001;

L=10; %length of the domain
nodesx = round(L/dx)+1.00; %number of nodes

%% Initialise k

k=zeros(nodesx,1);

%k heterogeneous - quadratic u shaped k(x)=1+ 0.1*(x-5)^2
k_init_condition_function=@(x) 1+ 0.1*(x-5)^2;
k_is_homogeneous=0;

for i=1:nodesx
    k(i) = k_init_condition_function((i-1)*dx);
end

%% Initialise a

a=zeros(nodesx,1);

%a heterogeneous lin increasing a, gradient 0.05
a_init_condition_function=@(x)  0.05*x;
a_is_homogeneous=0;

for i=1:nodesx
    a(i) = a_init_condition_function((i-1)*dx);
end

%% Eta

eta=1;

%% Initial density

init_condition_function = @(x) (gaussian_n/(sqrt(2*pi*gaussian_sigma^2)))*(exp(-((x-gaussian_mu).^2)/(2*gaussian_sigma^2)));

%Density
q0 = zeros(nodesx,1);
for i=1:nodesx
    q0(i) = init_condition_function((i-1)*dx);
end


%% Initial time and timestop

t0=0;
timestop=max(time_vector_record)+0.5;

%% Store initial variables (and initialise matrix)
max_recorded_data_points = 100000; %value to initialise for below
loop_count_stored = 1; %counter for when timestep data is recorded

q_hist = zeros(nodesx,max_recorded_data_points);
q_hist(:,loop_count_stored) = q0;

k_hist = zeros(nodesx, max_recorded_data_points);
k_hist(:,loop_count_stored) = k;

a_hist = zeros(nodesx, max_recorded_data_points);
a_hist(:,loop_count_stored) = a;

node_velocity_hist=zeros(nodesx,max_recorded_data_points);

dt_hist = zeros(nodesx, max_recorded_data_points);

t_hist= zeros(nodesx, max_recorded_data_points);

%% Node velocity

node_velocity = zeros(nodesx,1);
node_velocity_sign = zeros(nodesx,1);
for i=2:(nodesx-1)
    node_velocity(i) = (1/q0(i))*(1/eta)*(     k(i+1)*((1/(q0(i+1)))-a(i+1) ) -  k(i-1)*((1/(q0(i-1)))-a(i-1) )  );
    node_velocity_sign(i) = sign(node_velocity(i));
end

% Store the variables to use in the temporal loop
q_prev=q0;


%% When to record

min_diff_to_time_record_value = 0.001;

%% Time loop

t=t0;
loop_count=0;


while t < timestop
    
    loop_count=loop_count+1;
    
    %% Adaptive Timestepping
    min_diff_to_time_record = min(abs(time_vector_record -t));
    dt_new = (dt_courant*dx^2)/max(abs(node_velocity));
    if (min_diff_to_time_record < min_diff_to_time_record_value)
        dt =   min([dt_near_time_record,dt_ori*dx^2,10^(round(log10(dt_new),0))]);
    else
        dt =   min(dt_ori*dx^2,10^(round(log10(dt_new),0)));
    end
    dt_dxdx = dt/(dx^2);
    
    %% Update time
    t = t + dt;
    
    %% Solving for q
    
    %Can solve for q(1) using q(2) formula
    %Can solve for q(nodesx) using q(nodesx - 1) formula
    %Want to create a tridiagonal matrix to solve for q(2), ..., q(nodesv - 2),q(nodesv - 1).
    
    %Create the matrix to be used for the explicit difference iteration
    %subdiagonal
    subdiagonalq = zeros(nodesx,1);
    for i=2:nodesx
        if i==(nodesx-1)
            subdiagonalq(i) = -dt_dxdx*0.5*(k(i-1)/eta)*(1/q_prev(i-1)^2)*(2/3);
        else
            subdiagonalq(i) = -dt_dxdx*0.5*(k(i-1)/eta)*(1/q_prev(i-1)^2);
        end
    end
    
    %diagonal
    diagonalq = zeros(nodesx,1);
    for i=1:nodesx
        if i==2
            diagonalq(i) = -1 + dt_dxdx*0.5*(k(i)/eta)*(1/q_prev(i)^2)*(2/3);
        elseif i==(nodesx -1 )
            diagonalq(i) = -1 + dt_dxdx*0.5*(k(i)/eta)*(1/q_prev(i)^2)*(2/3);
        else
            diagonalq(i) = -1 + 2*dt_dxdx*0.5*(k(i)/eta)*(1/q_prev(i)^2);
        end
    end
    
    %superdiagonal (the first (nodesu -1) entries are the superdiagonal entries for nodes 1, 2,..., (nodesu-1)
    superdiagonalq = zeros(nodesx,1);
    for i=1:nodesx-1     %for spdiags
        if i==2
            superdiagonalq(i) = -dt_dxdx*0.5*((k(i+1)/eta)*(1/q_prev(i+1)^2)*(2/3));
        else
            superdiagonalq(i) = -dt_dxdx*0.5*(k(i+1)/eta)*((1/q_prev(i+1))^2);
        end
    end
    
    %Create the constant vector due to the resting spring length a
    constant_a_matrix_q = zeros(nodesx,1);
    for i=2:nodesx-1
        if i==2
            constant_a_matrix_q(i) = 0.5*dt_dxdx*( ((k(i+1)/eta)*(-a(i+1))) -2*((k(i)/eta)*(-a(i)))  )  + 0.5*dt_dxdx*(  -(1/3)*(k(i+1)/eta)*(-a(i+1)) + (4/3)*(k(i)/eta)*(-a(i)))  + 0.5*dt_dxdx*( ((k(i+1)/eta)*((1/q_prev(i+1))-a(i+1))) -2*((k(i)/eta)*((1/q_prev(i))-a(i))) + ((k(i-1)/eta)*((1/q_prev(i-1))-a(i-1))));
        elseif i == nodesx - 1
            constant_a_matrix_q(i) = 0.5*dt_dxdx*( -2*((k(i)/eta)*(-a(i))) + ((k(i-1)/eta)*(-a(i-1)))) +  0.5*dt_dxdx*(   (4/3)*(k(i)/eta)*(-a(i)) - (1/3)*(k(i-1)/eta)*(-a(i-1)))   + 0.5*dt_dxdx*( ((k(i+1)/eta)*((1/q_prev(i+1))-a(i+1))) -2*((k(i)/eta)*((1/q_prev(i))-a(i))) + ((k(i-1)/eta)*((1/q_prev(i-1))-a(i-1)))) ;
        else
            constant_a_matrix_q(i) = 0.5*dt_dxdx*( ((k(i+1)/eta)*(-a(i+1))) -2*((k(i)/eta)*(-a(i))) + ((k(i-1)/eta)*(-a(i-1))))  + 0.5*dt_dxdx*( ((k(i+1)/eta)*((1/q_prev(i+1))-a(i+1))) -2*((k(i)/eta)*((1/q_prev(i))-a(i))) + ((k(i-1)/eta)*((1/q_prev(i-1))-a(i-1)))) ;
        end
    end
    
    %Using function_06_tridia2 function (Thomas algorithm)
    implicit_resultq = function_06_tridia2(subdiagonalq(3:nodesx-1), diagonalq(2:nodesx-1), superdiagonalq(2:nodesx-2), (-q_prev(2:end-1) + constant_a_matrix_q(2:end-1) ));
    
    %includes the formulas used to determine the value of q at node 1 and the value of q at node nodesx
    %Result using the implicit method
    q = [1/( a(1) + (eta/k(1))*((1/3)*( (4*k(2)/eta)*( (1/implicit_resultq(1)) -a(2)) - (k(3)/eta)*( (1/implicit_resultq(2)) -a(3))     ))  )   ;
        implicit_resultq(1:end);
        1/( a(nodesx) + (eta/k(nodesx))*((1/3)*( (4*k(nodesx -1)/eta)*( (1/implicit_resultq(end)) -a(nodesx-1)) - (k(nodesx -2 )/eta)*( (1/implicit_resultq(end -1)) -a(nodesx-2))     ))  )  ];
    
    %% Update the node velocity
    
    %Determine the value of the velocity at each internal node 2,3,..., nodesv -1
    node_velocity = zeros(nodesx,1);
    node_velocity_sign = zeros(nodesx,1);
    for i=2:(nodesx-1)
        node_velocity(i) = (1/q(i))*(1/eta)*(     (k(i+1)*((1/(q(i+1)))-a(i+1) )) -  (k(i-1)*((1/(q(i-1)))-a(i-1)) )      );
        node_velocity_sign(i) = sign(node_velocity(i));
    end
    
    %% Update the spring stiffness function using the evolution equation for k
    
    %First store the previous value of k
    k_prev = k;
    if k_is_homogeneous ==1
    else
        %Now update using the new value of k
        k(1) = k_prev(1);
        k(nodesx) = k_prev(nodesx);
        for j=2:(nodesx - 1)
            if node_velocity_sign(j) == 1      %positive velocity - use information from the left (backward difference)
                k(j) = k_prev(j) - dt_dxdx*(  (1/q(j))*(  ((k(j)/eta)*((1/q(j))-a(j))) - ((k(j-1)/eta)*((1/q(j-1))-a(j-1)))  )*(k(j) -k(j-1))) ;
            elseif node_velocity_sign(j) == -1 %negative velocity - use infromation from the right (forward difference)
                
                k(j) = k_prev(j) - dt_dxdx*( (1/q(j))*(  ((k(j+1)/eta)*((1/q(j+1))-a(j+1))) - ((k(j)/eta)*((1/q(j))-a(j))) )*(k(j+1) -k(j)));
            else
                k(j) = k_prev(j) ;
            end
        end
        
    end
    
    
    %% Update the resting cell length function using the evolution equation for a
    
    a_prev=a;
    if a_is_homogeneous == 1
    else
        a(1) = a_prev(1);
        a(nodesx) = a_prev(nodesx);
        for j=2:(nodesx - 1)
            if node_velocity_sign(j) == 1      %positive velocity - use information from the left (backward difference)
                a(j) = a_prev(j) - dt_dxdx*(  (1/q(j))*(  ((k(j)/eta)*((1/q(j))-a(j))) - ((k(j-1)/eta)*((1/q(j-1))-a(j-1)))  )*(a(j) -a(j-1))) ;
            elseif node_velocity_sign(j) == -1 %negative velocity - use infromation from the right (forward difference)
                a(j) = a_prev(j) - dt_dxdx*( (1/q(j))*(  ((k(j+1)/eta)*((1/q(j+1))-a(j+1))) - ((k(j)/eta)*((1/q(j))-a(j)))  )*(a(j+1) -a(j)));
            else
                a(j) = a_prev(j) ;
            end
        end
        
    end
    
    %% Update q_prev
    
    q_prev =q;
    
    %% store q if want to store the value
    if  (min_diff_to_time_record < min_diff_to_time_record_value)
        loop_count_stored = loop_count_stored + 1;
        q_hist(:,loop_count_stored) = q;
        k_hist(:,loop_count_stored) = k;
        a_hist(:,loop_count_stored) = a;
        node_velocity_hist(:,loop_count_stored) = node_velocity;
        dt_hist(loop_count_stored) = dt;
        t_hist(loop_count_stored) = t;
    end
    t
    
end

%% Save
save([filepath_save '\PDE_results.mat'],'-v7.3')

end

