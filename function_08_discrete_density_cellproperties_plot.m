function function_08_discrete_density_cellproperties_plot(k, a,L, eta_cell,filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)
%function_discrete_density_plot

%k - cell stiffness for each cell
%a - resting cell length for each cell
%L - length of the domain
%filepath_save_figs - filepath to save figures
%time_vector_dots - time to plot lines with dots
%soln_discrete_Nm -
%N_cells - number of cells
%m_springs_per_cell - number of springs per cell


total_springs = N_cells*m_springs_per_cell;

if ismember(colouring,[1,2,3]) == 1 %for density, cell stiffness or resting cell length
    
    plot_resolution = 10;
    x_plot = zeros((total_springs)*plot_resolution,1);
    y_plot = zeros((total_springs)*plot_resolution,1);
    
    count=0;
    figure
    y=transpose(1:1:total_springs+1);
    
    for i=time_vector_dots
        count =count  +1;
        x=deval(soln_discrete_Nm, i);
        
        if colouring == 1
            dy=(1/m_springs_per_cell)*diff(y)./diff(x); % density
        elseif colouring == 2
            % cell stiffness
            dy=k;
        elseif colouring ==3
            % resting cell length
            dy=a;
        end
        
        
        %dy is the value in the intervals convert this so that
        %it can be plotted.
        
        for j=1:1:(total_springs)
            %New Positions dependent on the resolution
            
            difference = abs(x(j+1) - x(j));
            
            for r=1:1:plot_resolution
                x_plot((j-1)*plot_resolution + r) = x(j) + difference*(r-1)/(plot_resolution-1);
            end
            
            if j==1
                for r=1:1:plot_resolution
                    y_plot((j-1)*plot_resolution + r) = dy(j);
                end
            elseif j==total_springs
                for r=1:1:plot_resolution
                    y_plot((j-1)*plot_resolution + r) = dy(j);
                end
                
            else
                for r=1:1:plot_resolution
                    y_plot((j-1)*plot_resolution + r) = dy(j);
                end
            end
        end
        
        stairs(x_plot,y_plot);
        
        hold on
        legendInfo{count} = ['t = ' num2str(i)];
        
        xlabel('x')
        legend(legendInfo, 'location', 'south','Orientation','vertical')
        xticks([0,2.5,5,7.5,10])
        
        
        if colouring == 1
            %density
            ylabel('q')
            print(gcf,'-depsc2',[filepath_save_figs '\' 'Density_Multiple' '_' num2str(total_springs) '.eps']);
            saveas(gcf, [filepath_save_figs '\' 'Density_Multiple' '_' num2str(total_springs) '.fig'])
        elseif colouring == 2
            %cell stiffness
            ylabel('k')
            print(gcf,'-depsc2',[filepath_save_figs '\' 'Cellstiffness_Multiple_' num2str(total_springs) '.eps']);
            saveas(gcf, [filepath_save_figs '\' 'Cellstiffness_Multiple_' num2str(total_springs) '.fig'])
        elseif colouring == 3
            %resting cell length
            ylabel('a')
            print(gcf,'-depsc2',[filepath_save_figs '\' 'Restingcelllength_Multiple_' num2str(total_springs) '.eps']);
            saveas(gcf, [filepath_save_figs '\' 'Restingcelllength_Multiple_' num2str(total_springs) '.fig'])
        end
        
        
        
        
        
        
    end
    
elseif colouring == 4 % velocity
    
    %if greater than max v
    time_vector_dots(time_vector_dots - max(soln_discrete_Nm.x) >= 0) = max(soln_discrete_Nm.x);
    
    %determine the velocities of the spring boundaries
    velocity_u = zeros(total_springs+1,size(time_vector_dots,2));
    
    for j=1:1:size(time_vector_dots,2)
        for i=1:1:total_springs+1
            x=deval(soln_discrete_Nm,time_vector_dots(j));
            %first particle
            if i==1
                velocity_u(i,j) = 0; % zero velocity at fixed end
            elseif i==total_springs+1%end particles
                velocity_u(i,j) = 0; % zero velocity at fixed end
                
            else %velocity of spring boundaries
                velocity_u(i,j) = (1/(eta_cell/m_springs_per_cell))*(-(k(i-1)*m_springs_per_cell)*(x(i)-x(i-1) -a(i-1)/m_springs_per_cell) + (k(i)*m_springs_per_cell)*(x(i+1)-x(i) -a(i)/m_springs_per_cell));
                
            end
        end
    end
    
    %plot the figure
    
    figure
    for j=1:1:size(time_vector_dots,2)
        x=deval(soln_discrete_Nm,time_vector_dots(j));
        plot(x,velocity_u(:,j),'--')
        legendInfo{j} = ['t = ' num2str(time_vector_dots(j))];
        hold on
    end
    legend(legendInfo,'location', 'south')
    xlabel('x')
    ylabel('u')
    xticks([0,2.5,5,7.5,10])
    
    
    %save the figures
    print(gcf,'-depsc2',[filepath_save_figs '\' 'Velocity_Multiple_' num2str(total_springs) '.eps']);
    saveas(gcf, [filepath_save_figs '\' 'Velocity_Multiple_' num2str(total_springs) '.fig'])
    
end

