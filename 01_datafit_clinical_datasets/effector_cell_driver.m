%-----------------------------------------------------------------------------
% xCellignce data fitting - effector cells (CAR T-cells) and target cells
% (tumor cells)
% Created by: Viren S. and Justin W.
% Last update date: 1/8/2025
% Prompt: effector cell data fit
%-----------------------------------------------------------------------------%

% This driver function is meant to run the monte carlo simulation after
% preprocessing has been done.
% Before running run "effector_cell_data.m".
% All functions used are located in "param_est.m".
% As written, this script will generate a "results" folder with plots of
% the top 5 fits and a csv off all the fits with their loss.

run("effector_cell_data.m")

%% Set Simulation Options
warning('off', 'all');

%% Set model parameters --------------------------------------------------------

% Effector - m1 (values here no longer relevant for uniform monte carlo run)
K = 5e5;	        %Effector cell carrying capacity (cells)
k = 1e4;	        %Steepness of effector cell recruitment curve (cells*hours\(^{-1}\))
kc = 0.10;	        %Saturation level of tumor kill by CAR T-cells (hours\(^{-1}\))
kp2 = 0.01;	        %Max proliferation rate of CAR T-cells by tumor lysis (hours\(^{-1}\))
kd = e_kd;	        %CAR T-cell death rate (hour\(^{-1}\))
l = 1.00;           %ratio dependence (unitless)
s = 0.01;           %steepness (unitless)

pars = [K, k, kc, kp2, kd, l, s];
p_names = ["K", "k", "kc", "kp2", "kd", "l", "s"];

%Parameter lower/upper bounds - m1 - open
lb = [6.5e4 1e2 0.01 0.01 0.0001 0.01 0.01]; 
ub = [7.8e5 4e4 0.35 0.15 0.5000 3.50 3.50]; 

% Initial Concentrations (Assumes target_ic from previous run file)
c0 = [target_ic 0 0 0 0;
      target_ic 0 0 0 0;
      target_ic 0 0 0 0];

% Initial concentration of effector cells for three ratios
effector0 = [((25.0 * target_ic) * te)
             ((6.25 * target_ic) * te)
             ((1.00 * target_ic) * te)];

% Initial concentration of T cells for three ratios
tc0 = [((25.0 * target_ic) * (1-te))
       ((6.25 * target_ic) * (1-te))
       ((1.00 * target_ic) * (1-te))];

% Create matrix of parameter estimation data (from previous run file)
mean_targets = [mean_e25t1 mean_e6t1 mean_e1t1];
std_targets = [std_e25t1 std_e6t1 std_e1t1];
slope = m; intercept = b; p_log = parameters_logfit; tspan = time_data;

%% Perform simulations and save data.
% Perform Monte Carlo optimization

[all_pars, all_losses, init_pars, best_par] = param_est.monte_carlo_optimization(...
    num_simulations, pars, lb, ub, c0, tspan, mean_targets, ...
    p_log, effector0, tc0, slope, intercept, qfit, Ein, et, sampling_method, ...
    use_parallel, num_grid_points);
disp(all_losses)





