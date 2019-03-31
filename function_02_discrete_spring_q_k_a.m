function function_02_discrete_spring_q_k_a(  filepath_save )

%discrete_spring_q_k_a - run the discrete simulation for specific initial conditions for q,k,a and varying m_springs_per_cell

%% Define the initial conditions and setttings
reltol=1.0e-10;
abstol=1.0e-10;
odeeeoptions=odeset('RelTol',reltol,'AbsTol',abstol,'Events',@function_04_eventfun);

global ncall

global gaussian_n
global gaussian_sigma
global gaussian_mu
global eta_spring
global m_springs_per_cell
global num_springs
global L
global k_spring
global a_spring

%Initial conditions for the density condition to be a normal distribution
gaussian_n=11.056818502502322;
gaussian_sigma=sqrt(9);
gaussian_mu=5;


for m_springs_per_cell=[1,2,4]
    
    timestop=5000; %time to finish running the discrete simulation
  
    N_cells=10; %number of cells
    num_springs=N_cells*m_springs_per_cell; %number of springs
    L=10; %length of the domain
    eta_cell=1; %viscosity coefficient for a cell
    eta_spring = eta_cell/m_springs_per_cell; %viscosity
    
    global xl
    global xu
    xl=0;
    xu=L;
        
    %% Convert initial density into positions
    
    x0 = linspace(0,L,num_springs+1);
    x= fsolve(@function_03_fsolve_initial_positions_gaussian,x0)';
    pos=x';
    
    %% Convert initial spring stiffness and resting spring length into spring properties
    k_spring=zeros(num_springs,1);
    
    %k heterogeneous - quadratic u shaped k(x)=1+ 0.1*(x-5)^2
    function_init_k_spring_stiffness = @(x) 1+ 0.1*(x-5)^2;
    
    for i=1:num_springs
        k_posa =pos(i);
        k_posb =pos(i+1);
        k_gradm = 0.1;
        k_yconstant = 1;
        k_x0= [pos(i), pos(i+1)];
        k_quad_centre = 5;
        
        k_fun_k_stiffness_1 = @(x,k_posa,k_posb,k_gradm,k_yconstant,k_quad_centre) 2*( k_gradm*((x-k_quad_centre).^3)/3 + k_yconstant*x) - (k_gradm*((k_posa-k_quad_centre).^3)/3 + k_yconstant*k_posa) - ( k_gradm*((k_posb-k_quad_centre).^3)/3 + k_yconstant*k_posb);
        k_fun_x_k_stiffness_1  = @(x) k_fun_k_stiffness_1(x,k_posa,k_posb,k_gradm,k_yconstant,k_quad_centre);
        k_pos_to_plot(i)  = fzero(k_fun_x_k_stiffness_1,k_x0);
    end
    
    for i_loop2=1:num_springs
        k_spring(i_loop2) = function_init_k_spring_stiffness(k_pos_to_plot(i_loop2) )*(m_springs_per_cell);
    end
    
    
    
    a_spring=zeros(num_springs,1);
    
    %a heterogeneous lin increasing a, gradient 0.05
    function_init_a_resting_spring = @(x) 0.05*x;
    
    for i=1:num_springs
        a_posa =pos(i);
        a_posb =pos(i+1);
        a_gradm = 0.05;
        a_yintercept = 0;
        a_x0= [pos(i), pos(i+1)];
        a_fun_k_stiffness_1 = @(x,a_posa,a_posb,a_gradm,a_yintercept) 2*( (a_gradm*x.^2)/2 + a_yintercept*x) - ((a_gradm*a_posa^2)/2 + a_yintercept*a_posa) - ( (a_gradm*a_posb^2)/2 + a_yintercept*a_posb);
        a_fun_x_k_stiffness_1  = @(x) a_fun_k_stiffness_1(x,a_posa,a_posb,a_gradm,a_yintercept);
        a_pos_to_plot(i)  = fzero(a_fun_x_k_stiffness_1,a_x0);
    end
    
    for i_loop2=1:num_springs
        a_spring(i_loop2) = function_init_a_resting_spring(a_pos_to_plot(i_loop2) )/(m_springs_per_cell);
    end
    
    %% Use function to obtain the general solution
    
    % Independent variable for ODE integration
    ncall=0;
    
    soln_discrete = ode15s(@function_02_01_discrete_odes_for_ode15s,[0,timestop],pos,odeeeoptions);
    time_end =  max(soln_discrete.x);
    
    %Store the results for the different numbers of springs per cell
    %rescaled back to the cell level problem
    if m_springs_per_cell==1
        soln_discrete_m1 = soln_discrete;
        k_m1 = k_spring/(m_springs_per_cell);
        a_m1 = a_spring*(m_springs_per_cell);
        time_end_m1 = time_end;
    elseif m_springs_per_cell==2
        soln_discrete_m2 = soln_discrete;
        k_m2 = k_spring/(m_springs_per_cell);
        a_m2 = a_spring*(m_springs_per_cell);
        time_end_m2 = time_end;
    elseif m_springs_per_cell==4
        soln_discrete_m4 = soln_discrete;
        k_m4 = k_spring/(m_springs_per_cell);
        a_m4 = a_spring*(m_springs_per_cell);
        time_end_m4 = time_end;
    end
    
end

save([filepath_save '\DISCRETE_ALL_VARIABLES.mat'],'-v7.3')

end