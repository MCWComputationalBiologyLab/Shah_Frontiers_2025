
classdef param_est
    methods(Static)

%% Monte Carlo Optimization

function [all_pars, all_losses, init_pars, best_par] = monte_carlo_optimization(...
    num_simulations, p, lb, ub, u0, tspan, data, ...
    p_log, effector0, tc0, slope, intercept, q, Ein, et, ...
    sampling_method, use_parallel, num_grid_points)

% Run Monte Carlo parameter estimation or grid search with an option to run in parallel.
%
% Inputs: 
%   num_simulations - Number of simulations to run
%   p               - Parameter vector (used as mean for Gaussian sampling)
%   lb              - Lower boundary parameter vector
%   ub              - Upper boundary parameter vector
%   u0              - Initial concentration vector
%   tspan           - Time to run model 
%   data            - Comparison data for parameter estimation
%   p_log           - Log fit parameters (a,b)
%   effector0       - Initial concentrations for effector data
%   tc0             - Initial concentrations for t cell data
%   slope           - Slope value for conversion from CI to cell
%   intercept       - Intercept value for conversion of CI to cell
%   q               - Fit value used in cell calculation
%   Ein             - Index of effector cells added to experiment 
%   et              - Length of time data 
%   sampling_method - String to specify 'uniform', 'gaussian', or 'grid'
%   use_parallel    - Boolean flag to indicate whether to run in parallel
%   num_grid_points - Number of grid points for grid search (optional, default: 15)

% Outputs:
%   all_pars  - All parameter vectors estimated
%   all_losses - All corresponding loss values for all simulations
%   init_pars - Parameter vectors initialized using random vars 
%   best_par  - Most minimized set of parameters

    % Set default number of grid points if not provided
    if nargin < 1
        num_grid_points = 15;  % Default value for grid points
    end

    % Initialize arrays to store results
    all_pars = zeros(num_simulations, length(p));
    init_pars = zeros(num_simulations, length(p));
    all_losses = zeros(num_simulations, 1);
    all_combinations = 0;

    % Handle grid search setup if applicable
    if strcmp(sampling_method, 'grid')
        % Generate grid points for each parameter
        grid = cell(length(lb), 1);
        for i = 1:length(lb)
            grid{i} = linspace(lb(i), ub(i), num_grid_points);
        end
        
        % Generate all possible combinations of grid points
        [grid_combinations{1:length(grid)}] = ndgrid(grid{:});
        grid_combinations = cellfun(@(x) x(:), grid_combinations, 'UniformOutput', false);
        all_combinations = [grid_combinations{:}]; % All possible parameter sets
        num_simulations = size(all_combinations, 1) % Update to total number of grid points
        disp(num_simulations)
    end
    
    % Monte Carlo or grid search loop
    if use_parallel
        parfor i = 1:num_simulations
            warning('off', 'all')
            % Display progress
            disp(['Running simulation ', num2str(i), ' out of ', ...
                num2str(num_simulations)]);
            
            % Choose the appropriate sampling method
            initial_guess = param_est.select_sampling_method(...
                sampling_method, lb, ub, p, i, all_combinations);
            
            % Run the optimization
            [best_par, best_loss] = param_est.optimize_once(...
                initial_guess, lb, ub, tspan, u0, data, ...
                p_log, effector0, tc0, slope, intercept, q, Ein, et);
        
            % Store the results
            all_pars(i, :) = best_par;
            all_losses(i) = best_loss;
            init_pars(i, :) = initial_guess;
            disp(['Simulation ', num2str(i),' complete']);
        end
    else
        for i = 1:num_simulations
            % Display progress
            disp(['Running simulation ', num2str(i), ' out of ',...
                num2str(num_simulations)]);
            
            % Choose the appropriate sampling method
            initial_guess = param_est.select_sampling_method(...
                sampling_method, lb, ub, p, i, all_combinations);
            
            % Run the optimization
            [best_par, best_loss] = param_est.optimize_once(...
                initial_guess, lb, ub, tspan, u0, data, ...
                p_log, effector0, tc0, slope, intercept, q, Ein, et);
        
            % Store the results
            all_pars(i, :) = best_par;
            all_losses(i) = best_loss;
            init_pars(i, :) = initial_guess;
            disp(['Simulation ', num2str(i),' complete']);
        end
    end

    % Find the best parameters
    [~, best_idx] = min(all_losses);
    best_par = all_pars(best_idx, :);

end

%% Select sampling method

function initial_guess = select_sampling_method(sampling_method, lb, ub, p, ...
        i, all_combinations)

    % Helper function to select the initial guess based on sampling method
    switch sampling_method
        case 'uniform'
            % Standard Monte Carlo: random initial guess from uniform distribution
            initial_guess = lb + rand(size(lb)) .* (ub - lb);

        case 'gaussian'
            % Monte Carlo: random initial guess from Gaussian distribution 
            % with 'p' as mean
            initial_guess = normrnd(p, (ub - lb) / 6); 

            % Ensure the guess stays within bounds
            initial_guess = max(min(initial_guess, ub), lb);

        case 'grid'
            % Select initial guess from the pre-generated grid combinations
            initial_guess = all_combinations(i, :);   

    end
end

%% Estimation function for a single initial guess

function [par, loss] = optimize_once(p, lb, ub, tspan, c0, mean_targets, ...
    p_log, effector0, tc0, slope, intercept, q, Ein, et)
    
    % This function is meant to run fmincon once on a set of parameters.
    % Is iterated each Monte Carlo Simulation. 

    % Set Options for fmincon
    options = optimoptions('fmincon', ...
    'TolFun', 1e-12, ...
    'TolX', 1e-12, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'UseParallel','always', ...
    'Display', 'off');       % To suppress fmincon output

    % fmincon trycatch

    try
        % Perform Optimization using fmincon
        p = log(p); lb = log(lb); ub = log(ub); % conversion to log space

        [par, loss] = fmincon(@(x) param_est.lossFunction(x, tspan, c0, mean_targets, ...
            p_log, effector0, tc0, slope, intercept, q, Ein, et), ...
            p, [], [], [], [], lb, ub, [], options);

        par = exp(par);
        
    catch ME
        % If there's an error in fmincon or the model, return Inf as the loss
        disp('Optimization failed due to:');
        disp(ME.message);
        par = p; % Return the original parameter guess
        loss = Inf; % Assign a very high loss value so this run is not selected
    end

end

%% SSR function
      
function loss = lossFunction( p, tspan, c0, mean_CI, p_log, effector0, tc0, slope,...
        intercept, q, Ein, et)

    % Computes the sum of squared residuals between model and experimental data.
    %
    %  INPUTS:
    %     p         - (Array) Parameters to be estimated in the model.
    %     tspan     - (Array) Time vector used for ODE simulations.
    %     c0        - (Matrix) Initial concentrations for the target data
    %                 corresponding to different effector trials.
    %     mean_CI   - (Matrix) Mean CI data from three effector trials (e25, e6, e1).
    %     p_log     - (Array) Log fit parameters (used within the model function).
    %     effector0 - (Array) Initial effector concentrations for the simulation.
    %     slope     - (Scalar) Slope for the linear conversion from CI to cells.
    %     intercept - (Scalar) Intercept for the conversion from CI to cells.
    %     q         - (Scalar) Fit parameter used in CI to cell conversion.
    %     Ein       - (Integer) Index indicating the time when effector cells
    %                 are added to the experiment.
    %     et        - (Integer) Timepoint marking the end of the experimental data.
    %
    %   OUTPUT:
    %     loss - (Scalar) Sum of squared residuals (SSR) between model-predicted
    %            and experimental CI data, which is to be minimized during
    %            optimization.

    try
        % Set up the ODE model
        p = exp(p); % convert back to linear space
        model = @(t,y)(cyto(t, y, p, p_log)); % Define the model function
        t_tumor = tspan(1:Ein-1); % Time before effector addition
        t_effector = tspan(Ein:et); % Time after effector addition
    
        % Solve the ODEs for each effector trial (E25, E6.25, E1)
        [~, fit25_1] = ode15s(model, t_tumor, c0(1,:)');  % E25:T1
        [~, fit25_2] = ode15s(model, t_effector, [fit25_1(end,1) effector0(1) tc0(1) 0 0]');
    
        [~, fit6_1] = ode15s(model, t_tumor, c0(2,:)');   % E6.25:T1
        [~, fit6_2] = ode15s(model, t_effector, [fit6_1(end,1) effector0(2) tc0(2) 0 0]');
    
        [~, fit1_1] = ode15s(model, t_tumor, c0(3,:)');   % E1:T1
        [~, fit1_2] = ode15s(model, t_effector, [fit1_1(end,1) effector0(3) tc0(3) 0 0]');
    
        % Convert ODE solutions to cellular index (CI)
        CI_fit25 = param_est.convert_CI_to_cell(fit25_1, fit25_2, slope, intercept, Ein, q);
        CI_fit6  = param_est.convert_CI_to_cell(fit6_1, fit6_2, slope, intercept, Ein, q);
        CI_fit1  = param_est.convert_CI_to_cell(fit1_1, fit1_2, slope, intercept, Ein, q);
    
        % Calculate sum of squared residuals (SSR) between model and experimental data
        SSR_pts = ((mean_CI(:,1) - CI_fit25).^2) + ...
                  ((mean_CI(:,2) - CI_fit6).^2) + ...
                  ((mean_CI(:,3) - CI_fit1).^2);

        % Sum of all residuals (final loss value)
        loss = sum(SSR_pts, 'all');
    
    catch ME
        % Handle any errors during the ODE integration or other operations
        %disp("ODE solver failed. Moving to next iteration...");
        
        % Assign a large loss value in case of failure
        loss = Inf; % Large value to penalize unsuccessful runs

    end
end

%% Run model

function [fit_results, fit_results_tc, fit_results_ec, fit_results_utc, ...
          fit_results_fc, fit_results_fp] = runSingleParSet( ...
    p, ... % parameter estimate for effects
    tspan, ... % Time 
    c0, ... % Initial Concentrations for target data
    p_log, ... % parameter estimate for targets (a,b)
    effector0,...% Initial concentrations for effector data
    tc0,...% Initial concentrations for t cells
    slope, ... % Slope value for conversion for CI to cell
    intercept, ... % Intercept value for conversion of CI to cell
    q, ... % (???) some fit value used in cells calculation
    Ein, ... % Index of effector cells added to experiment 
    et ... % Length of time data (Can just calculate this in fxn but w/e)
    )
    
    try
        % Set Model
        model = @(t,y)(cyto(t, y, p, p_log));
        t_tumor = tspan(1:Ein-1);
        t_effector = tspan(Ein:et);
    
        % Solve for simulation data _________________________________________________
        % E25:T1
        [~, fit25_1] = ode15s(model, t_tumor, c0(1,:)');
        [~, fit25_2] = ode15s(model, t_effector, [fit25_1(end,1) effector0(1) tc0(1) 0 0]');

        % E6.25:T1 
        [~, fit6_1] = ode15s(model, t_tumor, c0(2,:)');
        [~, fit6_2] = ode15s(model, t_effector, [fit6_1(end,1) effector0(2) tc0(2) 0 0]');
    
        % E1:T1 
        [~, fit1_1] = ode15s(model, t_tumor, c0(3,:)');
        [~, fit1_2] = ode15s(model, t_effector, [fit1_1(end,1) effector0(3) tc0(3) 0 0]');
    
        % Perform Conversions
        CI_fit25 = param_est.convert_CI_to_cell(fit25_1, fit25_2, slope, intercept, Ein, q);
        CI_fit6  = param_est.convert_CI_to_cell(fit6_1, fit6_2, slope, intercept, Ein, q);
        CI_fit1  = param_est.convert_CI_to_cell(fit1_1, fit1_2, slope, intercept, Ein, q);

        fit_results = [CI_fit25, CI_fit6, CI_fit1];

        % Cells
        fit_results_tc(:,1) = [fit25_1(:,1); fit25_2(:,1)];
        fit_results_ec(:,1) = [fit25_1(:,2); fit25_2(:,2)];
        fit_results_utc(:,1) = [fit25_1(:,3); fit25_2(:,3)];
        fit_results_fc(:,1) = [fit25_1(:,4); fit25_2(:,4)];
        fit_results_fp(:,1) = [fit25_1(:,5); fit25_2(:,5)];
        fit_results_tc(:,2) = [fit6_1(:,1); fit6_2(:,1)];
        fit_results_ec(:,2) = [fit6_1(:,2); fit6_2(:,2)];        
        fit_results_utc(:,2) = [fit6_1(:,3); fit6_2(:,3)];
        fit_results_fc(:,2) = [fit6_1(:,4); fit6_2(:,4)];
        fit_results_fp(:,2) = [fit6_1(:,5); fit6_2(:,5)];
        fit_results_tc(:,3) = [fit1_1(:,1); fit1_2(:,1)];
        fit_results_ec(:,3) = [fit1_1(:,2); fit1_2(:,2)];
        fit_results_utc(:,3) = [fit1_1(:,3); fit1_2(:,3)];
        fit_results_fc(:,3) = [fit1_1(:,4); fit1_2(:,4)];
        fit_results_fp(:,3) = [fit1_1(:,5); fit1_2(:,5)];        
    
    catch ME
        % Handle any errors during the ODE integration or other operations
        disp("ODE solver failed. Moving to next iteration...")
        
        % Assign a large loss value in case of failure
        fit_results = 1e6; % You can adjust this value if needed
        
    end
end

%% CI to cell conversion
    
function CI_fit = convert_CI_to_cell(fit1, fit2, slope, intercept, Ein, q)

    % Extract fits
    tumor_fit = [fit1(:,1); fit2(:,1)];
    effector_fit = [fit1(:,2); fit2(:,2)] + [fit1(:,3); fit2(:,3)];

    
    % Convert to Cell
    tumor_fit_CI = (tumor_fit - intercept) ./ slope;
    effector_fit_CI(1:Ein-1) = 0; %There is no effector before this time
    effector_fit_CI(Ein:size(effector_fit,1)) = polyval(q,effector_fit(Ein:end));

    CI_fit = tumor_fit_CI + effector_fit_CI';

end

end
end
