%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Murphy2019-An individual-based mechanical model of cell movement in heterogeneous tissues and its coarse-grained approximation.
% 25/03/2019
% Code to generate Supplementary Figure 6 made available to GitHub
% Includes:
% - discrete model
% - pde model
% - characteristics figures with density, cell stiffness, resting cell length and velocity colouring
% - snapshot comparison (at one timepoint)
%       - density
%       - cell stiffness
%       - resting cell length
%       - velocity

%% Required functions to download include:

% function_02_discrete_spring_q_k_a.m
% function_03_fsolve_initial_positions_gaussian.m
% function_04_eventfun.m
% function_05_pde_spring_q_k_a.m
% function_06_tridia2.m
% function_07_characteristics_plot.m
% function_08_discrete_density_cellproperties_plot.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Model parameters

 %corresponds to normal density
 %heterogeneous k = 1+0.1(x-5)^2
 %heterogeneous a = 0.05x
filepath_save = 'C:\ResearchC\Paper 1\GitHub\Results';

%% Run the discrete model

function_02_discrete_spring_q_k_a(filepath_save)

%% Run the pde model
time_vector_record=[0,1.25,15];
function_05_pde_spring_q_k_a(  filepath_save, time_vector_record )

%% Load the discrete data

clear
filepath_load = 'C:\ResearchC\Paper 1\GitHub\Results';
load([filepath_load, '\DISCRETE_ALL_VARIABLES.mat']);

%% Figures - Characteristics parameters - for m_springs_per_cell=4

k=k_m4; %cell values 
a=a_m4;
time_start=0;
time_end=16.25;
time_vector_dots = [0,1.25,15];
soln_discrete_Nm = soln_discrete_m4;
N_cells=10;
m_springs_per_cell = 4;

xticks_vec=[0,2.5,5,7.5,10];
yticks_vec=[0,5,10,15];
plot_every_n_characteristics=1;
circle_size=35;

filepath_save_figs='C:\ResearchC\Paper 1\GitHub\Figs';

%% Figures - Characteristics - density colouring

colouring =1;
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Characteristics - cell stiffness colouring

colouring =2;
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Characteristics - resting cell length colouring

colouring =3;
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Characteristics - velocity colouring

colouring =4;
function_07_characteristics_plot(k, a,L, time_start, time_end, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring, eta_cell, xticks_vec, yticks_vec, plot_every_n_characteristics,circle_size)

%% Figures - Snapshots - density - discrete

colouring=1;
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Figures - Snapshots - cell stiffness - discrete

colouring=2;
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Figures - Snapshots - resting cell length - discrete

colouring =3;
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Figures - Snapshots - velocity - discrete

colouring =4;
function_08_discrete_density_cellproperties_plot(k, a,L,eta_cell, filepath_save_figs, time_vector_dots, soln_discrete_Nm, N_cells, m_springs_per_cell, colouring)

%% Load continuum data 

clear
filepath_load = 'C:\ResearchC\Paper 1\GitHub\Results';
load([filepath_load '\Results_PDE_results.mat']);

filepath_save_figs = 'C:\ResearchC\Paper 1\GitHub\Figs';

N_cell=10; %to open related figure from discrete results
m_springs_per_cell=4; %to open related figure from discrete results
eta_cell=eta;
time_vector_dots=[0,1.25,15];

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
