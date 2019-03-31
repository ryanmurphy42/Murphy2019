function function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)
%function_characteristcs_plot

%k - cell stiffness for each cell
%a - resting cell length for each cell
%L - length of the domain
%time_start - initial time to start characteristic plot
%time_end - final time for characteristic plot
%filepath_save_figs - filepath to save figures
%time_vector_dots - time to plot lines with dots
%soln_discrete_Nm -
%N_cells - number of cells
%m_springs_per_cell - number of springs per cell
%colouring - 1, Density, 2- Cell stiffness, 3-Resting cell length, 4-velocity
total_springs = N_cells*m_springs_per_cell;


set(0,'defaultAxesFontSize',18)

%determine times for the scatter plot
time_vector_characteristics = linspace(time_start,time_end,max(300,  max(1,round((time_end-time_start)/250)))  );

%if greater than max v
time_vector_characteristics(time_vector_characteristics - max(soln_discrete_Nm.x) >= 0) = max(soln_discrete_Nm.x);


%Add in an approximate particle density
particle_positions = zeros(total_springs+1,size(time_vector_characteristics,2));
for pp = time_vector_characteristics
    if    pp == min( time_vector_characteristics)
        particle_positions = deval(soln_discrete_Nm,pp);
    else
        particle_positions = [particle_positions, deval(soln_discrete_Nm,pp)];
    end
end

if colouring ==1
    %use the density as the colouring
    %Add in an approximate particle density for colouring
    particle_density = zeros(total_springs+1,size(time_vector_characteristics,2));
    loopcounter =0;
    for pp = time_vector_characteristics
        loopcounter = loopcounter+1;
        for i=1:total_springs+1
            if i ==total_springs+1
                particle_density(i,loopcounter) =  (1/m_springs_per_cell)*(1/(particle_positions(i,loopcounter) - particle_positions(i-1,loopcounter)));
            elseif i ==1
                particle_density(i,loopcounter) =  (1/m_springs_per_cell)*(1/(particle_positions(i+1,loopcounter) - particle_positions(i,loopcounter)));
            else
                particle_density(i,loopcounter) =   0.5*(1/m_springs_per_cell)*( (1/(particle_positions(i+1,loopcounter) - particle_positions(i,loopcounter))) + (1/(particle_positions(i,loopcounter) - particle_positions(i-1,loopcounter))));
            end
        end
    end
    
    % Add in a colour to show the density
    figure
    for m = 1:plot_every_n_characteristics:(total_springs+1)
        scatter(particle_positions(m,1:end),time_vector_characteristics(1:end),10, particle_density(m,1:end),'filled' )
        hold on
    end
    title('Characteristics - Density Colouring')
    
elseif colouring==2
    %use the cell stiffness for colouring
    particle_k = zeros(total_springs+1,size(time_vector_characteristics,2));
    loopcounter =0;
    for pp = time_vector_characteristics
        loopcounter = loopcounter+1;
        for i=1:total_springs+1
            if i ==total_springs+1
                particle_k(i,loopcounter) =  k(total_springs);
            elseif   i==1
                particle_k(i,loopcounter) =  k(1);
            else
                particle_k(i,loopcounter) = (k(i-1)+k(i))/2;
            end
        end
    end
    
    % Add in a colour to show the density
    figure
    for m = 1:plot_every_n_characteristics:(total_springs+1)
        scatter(particle_positions(m,1:end),time_vector_characteristics(1:end),10, particle_k(m,1:end),'filled' )
        hold on
    end
    title('Characteristics - Cell Stiffness Colouring')
    
elseif colouring ==3
    %use the resting cell length for colouring
    particle_a = zeros(total_springs+1,size(time_vector_characteristics,2));
    loopcounter =0;
    for pp = time_vector_characteristics
        loopcounter = loopcounter+1;
        for i=1:total_springs+1
            if i ==total_springs+1
                particle_a(i,loopcounter) =  a(total_springs);
            elseif   i==1
                particle_a(i,loopcounter) =  a(1);
            else
                particle_a(i,loopcounter) = (a(i-1)+a(i))/2;
            end
        end
    end
    
    % Add in a colour to show the density
    figure
    for m = 1:plot_every_n_characteristics:(total_springs+1)
        scatter(particle_positions(m,1:end),time_vector_characteristics(1:end),10, particle_a(m,1:end),'filled' )
        hold on
    end
    title('Characteristics - Resting Cell Length Colouring')
    
elseif colouring ==4
    %use the velocity for colouring
    Force_particles_at_time = zeros(total_springs+1,size(time_vector_characteristics,2));
    loopcounter =0;
    for pp = 1:size(time_vector_characteristics,2)
        loopcounter = loopcounter+1;
        for i=1:1:total_springs+1
            x=deval(soln_discrete_Nm,time_vector_characteristics(pp));
            if i==1
                Force_particles_at_time(i,loopcounter) = 0; %zero velocity imposed at boundary
            elseif i==total_springs+1
                Force_particles_at_time(i,loopcounter) = 0; %zero velocity imposed at boundary
            else %interior particles
                Force_particles_at_time(i,loopcounter) =(-(k(i-1)*m_springs_per_cell)*((x(i)-x(i-1)) -a(i-1)/m_springs_per_cell) + (k(i)*m_springs_per_cell)*((x(i+1)-x(i)) -a(i)/m_springs_per_cell));
            end
        end
    end
    
    particle_u = zeros(total_springs+1,size(time_vector_characteristics,2));
    loopcounter =0;
    for pp = time_vector_characteristics
        loopcounter = loopcounter+1;
        for i=1:total_springs+1
            if i ==total_springs+1
                %zero velocity at the boundary
                particle_u(i,loopcounter) =0;
            elseif i==1
                %zero velocity at the boundary
                particle_u(i,loopcounter) =0;
            else
                %(1/(eta))*(f_(i+1) - f_(i)) in the discrete model
                particle_u(i,loopcounter) = (1/(eta_cell/m_springs_per_cell))*(Force_particles_at_time(i,loopcounter)) ;
            end
        end
    end
    
    for m = 1:plot_every_n_characteristics:(total_springs+1)
        scatter(particle_positions(m,1:end),time_vector_characteristics(1:end),15, particle_u(m,1:end),'filled' )
        hold on
    end
    title('Characteristics - Velocity Colouring')
    
end


xlabel('x')
ylabel('t')
xlim([0,L])
xticks([0,2.5,5,7.5,10])
ylim([0,time_end])

% Add in lines for the snapshots
for i=1:size(time_vector_dots,2)
    x_plot_lines_cross = linspace(0,10,10);
    y_plot_lines_cross = time_vector_dots(i)*ones(10,1);
    plot(x_plot_lines_cross,y_plot_lines_cross,'linewidth',1.5,'color','k')
end

% Add in circle markers
for i=1:size(time_vector_dots,2)
    x_plot_lines_cross_springs = deval(soln_discrete_Nm,time_vector_dots(i));
    y_plot_lines_cross_springs = time_vector_dots(i)*ones(total_springs+1,1);
    
    x_plot_lines_cross_cells = x_plot_lines_cross_springs(1:plot_every_n_characteristics:total_springs+1);
    y_plot_lines_cross_cells = y_plot_lines_cross_springs(1:plot_every_n_characteristics:total_springs+1);
    
    scatter(x_plot_lines_cross_cells,y_plot_lines_cross_cells,circle_size,'filled','k')
end

%flip the y axis axis
ax = gca;
ax.YDir = 'reverse';

xticks(xticks_vec);
yticks(yticks_vec);


%save the figure with the colorbar
colorbar
print(gcf,'-depsc2',  [filepath_save_figs '\Characteristics_with_colourbar_' num2str(total_springs) '_type_' num2str(colouring) '.eps']);
saveas(gcf, [filepath_save_figs '\Characteristics_with_colourbar_' num2str(total_springs) '_type_' num2str(colouring) '.fig'])

end

