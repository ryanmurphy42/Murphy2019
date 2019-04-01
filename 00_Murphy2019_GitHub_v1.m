%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Murphy2019-An individual-based mechanical model of cell movement in heterogeneous tissues and its coarse-grained approximation.

% 31/03/2019

% Code to generate Supplementary Figure 6 made available to GitHub
% Includes:
% - running discrete model for 1,2,4 springs per cell
% - pde model
% - characteristics figures with density, cell stiffness, resting cell length and velocity colouring
% - snapshot comparison (at one timepoint)
%       - density
%       - cell stiffness
%       - resting cell length
%       - velocity

%% Required functions to have saved in matlab folder include:

% function_02_discrete_spring_q_k_a.m
% function_02_01_discrete_odes_for_ode15s.m
% function_03_fsolve_initial_positions_gaussian.m
% function_04_eventfun.m
% function_05_pde_spring_q_k_a.m
% function_06_tridia2.m
% function_07_characteristics_plot.m
% function_08_discrete_density_cellproperties_plot.m
% function_09_discrete_ctm_comparison_density_cellproperties_plot.m

%% Create 2 new folders if dont already exist

mkdir('Results') %folder to store the discrete and continuum simulation file
mkdir('Figs') %folder to save figures to

%% Reset matlab

clear
clc
close all

%% Model parameters

 %density - normal distribution - this is coded into the functions
 %heterogeneous k = 1+0.1(x-5)^2 - this is coded into the functions
 %heterogeneous a = 0.05x - this is coded into the functions
 
filepath_save = [pwd '\Results'];

%% Run the discrete model

function_02_discrete_spring_q_k_a(filepath_save)

%% Run the pde model

time_vector_record=[0,1.25,15]; %times at which want the continuum solution
function_05_pde_spring_q_k_a(  filepath_save, time_vector_record )

%% Load the discrete data

clear
filepath_load = [pwd '\Results'];
load([filepath_load, '\DISCRETE_ALL_VARIABLES.mat']);

%% Figures - Characteristics parameters - for m_springs_per_cell=4

k=k_m4; %stiffness at scaled back to cell level for 4 springs per cell
a=a_m4; %resting cell length scaled back to cell level for 4 springs per cell
time_start=0; %time to start the characteristic plot
time_end=16.25; %time to end the characteristic plot
time_vector_dots = [0,1.25,15]; %lines for times corresponding to snapshots
soln_discrete_Nm = soln_discrete_m4; %discrete solution result
N_cells=10; %number of cells in the system
m_springs_per_cell = 4; %number of springs per cell

xticks_vec=[0,2.5,5,7.5,10]; %x ticks (defined for L=10)
yticks_vec=[0,5,10,15]; %y ticks (chosen for the above time_start and time_end)
plot_every_n_characteristics=1; %number of characteristics plotted
circle_size=35; %size of dots to plot on lines

filepath_save_figs=[pwd '\Figs']; %location to save figures

%% Figures - Characteristics - density colouring

colouring =1; %corresponds to density
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Characteristics - cell stiffness colouring

colouring =2; %corresponds to cell stiffness
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Characteristics - resting cell length colouring

colouring =3; %correpsonds to resting cell length
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Characteristics - velocity colouring

colouring =4; %corresponds to velocity
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Snapshots - density - discrete

colouring=1; %corresponds to density
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Figures - Snapshots - cell stiffness - discrete

colouring=2; %corresponds to cell stiffness
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Figures - Snapshots - resting cell length - discrete

colouring =3; %corresponds to resting cell length
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Figures - Snapshots - velocity - discrete

colouring =4; %corresponds to velocity
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Load continuum data 

clear
filepath_load = [pwd '\Results'];
load([filepath_load '\PDE_results.mat']);

filepath_save_figs = [pwd '\Figs'];

N_cell=10; %to open related figure from discrete results
m_springs_per_cell=4; %to open related figure from discrete results
eta_cell=eta; %a relabelling
time_vector_dots=[0,1.25,15]; %times to plot snapshots

%% Figures - Snapshots - density - continuum comparison

colouring=1;
function_09_discrete_ctm_comparison_density_cellproperties_plot(t_hist, q_hist,k_hist, a_hist, eta_cell,L,dx, filepath_save_figs, time_vector_dots, N_cell, m_springs_per_cell,  colouring)

%% Figures - Snapshots - cell stiffness - continuum comparison

colouring=2;
function_09_discrete_ctm_comparison_density_cellproperties_plot(t_hist, q_hist,k_hist, a_hist, eta_cell,L,dx, filepath_save_figs, time_vector_dots, N_cell, m_springs_per_cell,  colouring)

%% Figures - Snapshots - resting cell length - continuum comparison

colouring =3;
function_09_discrete_ctm_comparison_density_cellproperties_plot(t_hist, q_hist,k_hist, a_hist, eta_cell,L,dx, filepath_save_figs, time_vector_dots, N_cell, m_springs_per_cell,  colouring)

%% Figures - Snapshots - velocity - continuum comparison

colouring =4;
function_09_discrete_ctm_comparison_density_cellproperties_plot(t_hist, q_hist,k_hist, a_hist, eta_cell,L,dx, filepath_save_figs, time_vector_dots, N_cell, m_springs_per_cell,  colouring)
