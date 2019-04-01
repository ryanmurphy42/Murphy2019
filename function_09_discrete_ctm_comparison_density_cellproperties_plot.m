function function_09_discrete_ctm_comparison_density_cellproperties_plot(t_hist, q_hist,k_hist, a_hist, eta_cell,L,dx, filepath_save_figs, time_vector_dots, N_cell, m_springs_per_cell,  colouring)
%function_09_discrete_ctm_comparison_density_cellproperties_plot - overlay the continuum results on top of the discrete results

total_springs = N_cell*m_springs_per_cell;

if colouring == 1 %density
    
    openfig([filepath_save_figs '\' 'Density_Multiple_' num2str(total_springs) '.fig']) %open the discrete result
    set(0,'DefaultLegendAutoUpdate','off') %legend does not need updating as it already states the 3 times
    
    %loop through times and plot the continuum result. Continuum solution must have been run for the same time points.
    tcomparing_count=0;
    for i=time_vector_dots
        tcomparing = i;
        tcomparing_count = tcomparing_count +1;
        [first_compare_diff first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        hold on
        plot(0:dx:L , q_hist(:,first_compare_index),'LineWidth',1, 'color','k')
        
    end
    
    %update figure properties
    ylabel('q')
    set(gca,'FontSize',18)
    
    %save figures
    saveas(gcf,[filepath_save_figs '\' 'Density_Multiple_Compared' num2str(total_springs) '.fig'])
    print(gcf,'-depsc2',[filepath_save_figs '\' 'Density_Multiple_Compared' num2str(total_springs) '.eps'])
    
elseif colouring  == 2 %cell stiffness
    
    openfig([filepath_save_figs '\' 'Cellstiffness_Multiple_' num2str(total_springs) '.fig'])
    set(0,'DefaultLegendAutoUpdate','off')
    
    tcomparing_count=0;
    for i=time_vector_dots
        tcomparing = i;
        tcomparing_count = tcomparing_count +1;
        [first_compare_diff first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        hold on
        plot(0:dx:L , k_hist(:,first_compare_index),'LineWidth',1, 'color','k')
        
    end
    
    ylabel('k')
    set(gca,'FontSize',18)
    
    saveas(gcf,[filepath_save_figs '\' 'Cellstiffness_Multiple_' num2str(total_springs) '.fig'])
    print(gcf,'-depsc2',[filepath_save_figs '\' 'Cellstiffness_Multiple_' num2str(total_springs) '.eps'])
    
elseif colouring == 3 %resting cell length
    
    openfig([filepath_save_figs '\' 'Restingcelllength_Multiple_' num2str(total_springs) '.fig'])
    set(0,'DefaultLegendAutoUpdate','off')
    
    tcomparing_count=0;
    for i=time_vector_dots
        tcomparing = i;
        tcomparing_count = tcomparing_count +1;
        [first_compare_diff first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        hold on
        plot(0:dx:L ,a_hist(:,first_compare_index),'LineWidth',1, 'color','k')
        
    end
    
    ylabel('a')
    set(gca,'FontSize',18)
    
    saveas(gcf,[filepath_save_figs '\' 'Restingcelllength_Multiple_' num2str(total_springs) '.fig'])
    print(gcf,'-depsc2',[filepath_save_figs '\' 'Restingcelllength_Multiple_' num2str(total_springs) '.eps'])
    
elseif colouring == 4 %velocity
    
    openfig([filepath_save_figs '\' 'Velocity_Multiple_' num2str(total_springs) '.fig'])
    set(0,'DefaultLegendAutoUpdate','off')
    
    tcomparing_count=0;
    for i=time_vector_dots
        tcomparing = i;
        tcomparing_count = tcomparing_count +1;
        [first_compare_diff first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        
        q_for_force_at_nodes_ctm_plot = q_hist(:,first_compare_index);
        k_for_force_at_nodes_ctm_plot = k_hist(:,first_compare_index);
        a_for_force_at_nodes_ctm_plot = a_hist(:,first_compare_index);
        
        force_at_nodes_ctm_plot = k_for_force_at_nodes_ctm_plot.*(1./(q_for_force_at_nodes_ctm_plot) -(a_for_force_at_nodes_ctm_plot));
        force_gradient_at_nodes_ctm_plot = (1./(eta_cell*q_for_force_at_nodes_ctm_plot)).*gradient( force_at_nodes_ctm_plot,dx);
        hold on
        plot(0:dx:L , force_gradient_at_nodes_ctm_plot(:),'LineWidth',1, 'color','k')
        
    end
    
    ylabel('u')
    set(gca,'FontSize',18)
    
    saveas(gcf,[filepath_save_figs '\' 'Velocity_Multiple_Compared' num2str(total_springs) '.fig'])
    print(gcf,'-depsc2',[filepath_save_figs '\' 'Velocity_Multiple_Compared' num2str(total_springs) '.eps'])
    
end


end

