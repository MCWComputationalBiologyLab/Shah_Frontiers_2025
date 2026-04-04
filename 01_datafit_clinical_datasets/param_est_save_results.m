
classdef param_est_save_results
    methods(Static)

%% Save results function

function [target_fit_CI,...
            sim_data, sim_data_tc, sim_data_ec,...
            sim_data_utc, sim_data_fc, sim_data_fp, R_squared, RMSE] ...
         = saveResults(...
         all_losses, all_pars, init_pars, tspan,...
         c0, p_log, effector0, tc0, ...
        mean_target_CI, std_target_CI, mean_target_c,...
        target_fit_CI, target_fit_log, ...
        exp_data, exp_data_std, columns,...
        cells_tar, CI_tar, cells_eff, CI_eff,...
        slope, intercept, q,...
        Ein, et, top_n, save_dir, p_names, csv_name, ...
        seed, method, makeResultFolder)
    % This function finds the top N parameter sets
    % with the lowest losses, runs the model with these parameters,
    % generates plots of the results, and creates a CSV file.
    %
    % Inputs:
    %   all_losses  - Array of loss values for all parameter sets.
    %   all_pars    - Matrix of parameter sets.
    %   init_pars   - Matrix of initial randomized parameter sets.
    %   tspan       - Time vector for the model.
    %   c0          - Initial concentrations vector.
    %   p_log       - Log fit parameters.
    %   effector0   - Initial concentrations for effector data.
    %   slope       - Slope value for conversion for CI to cell.
    %   intercept   - Intercept value for conversion of CI to cell.
    %   q           - Fit value used in cells calculation.
    %   Ein         - Index of effector cells added to experiment.
    %   et          - Length of time data.
    %   top_n       - Number of top results to consider.
    %   save_dir    - Directory to save the plots.
    %   p_names     - Cell array of parameter names for CSV headers.
    %   csv_name    - Name for the output CSV file.
    %   method      - Type of sampling used: uniform, gaussian, grid

    % Create the directory if it doesn't exist
    % if ~exist(save_dir, 'dir')
    %     mkdir(save_dir);
    % end

    % Ensure a unique result folder by appending a number if the folder already exists
    if makeResultFolder
        result_folder = fullfile(save_dir);
            counter = 1;
            while exist(result_folder, 'dir')
                result_folder = fullfile(save_dir, sprintf('%d', counter));
                counter = counter + 1;
            end
            mkdir(result_folder);
    else
        result_folder = fullfile(save_dir);
    end



    % Sort the all_losses array in ascending order
    [sorted_losses, sort_idx] = sort(all_losses);
    
    % Reorder all_pars to match the sorted all_losses
    sorted_pars = all_pars(sort_idx, :);
    
    % Select the top N parameters based on the lowest losses
    top_pars = sorted_pars(1:min(top_n, size(sorted_pars, 1)), :);
    
    % Loop through each set of top parameters
    for i = 1:size(top_pars, 1)
        try
            fprintf("Saving plot %d/%d\n", i, size(top_pars, 1))
            % Run the model for each top parameter set
            [sim_data, sim_data_tc, sim_data_ec,...
             sim_data_utc, sim_data_fc, sim_data_fp] ...
             = param_est.runSingleParSet(top_pars(i, :), tspan, c0,...
               p_log, effector0, tc0, slope, intercept, q, Ein, et);
            
            % Create a new figure for the current plot without displaying it
            fig = plots_indv_data(tspan, mean_target_CI, std_target_CI, ...
                                  mean_target_c, target_fit_CI, target_fit_log, ...
                                  exp_data, exp_data_std, sim_data, sim_data_tc, sim_data_ec, ...
                                  sim_data_utc, sim_data_fc, sim_data_fp, columns, Ein);

            % Save the plot in the result folder
            saveas(fig, fullfile(result_folder, sprintf('top_result_%d.png', i)));
            close(fig); % Close the figure to avoid cluttering memory
        catch ME
            fprintf('Error occurred while processing parameter set %d: %s\n', ...
                i, ME.message);
        end
    end

    % get to par data

    % Run the model for top parameter set
    [sim_data, sim_data_tc, sim_data_ec,...
     sim_data_utc, sim_data_fc, sim_data_fp] ...
     = param_est.runSingleParSet(top_pars(1, :), tspan, c0,...
       p_log, effector0, tc0, slope, intercept, q, Ein, et);

    % Calculate R^2
        
        % R^2 for 25:1

        y_bar = sum(exp_data(:,1) / length(exp_data(:,1)));
        SSY = sum((exp_data(:,1)  - y_bar).^2);
        SSE = sum((exp_data(:,1)-sim_data(:,1)).^2);
        R_squared(1) = round((SSY - SSE) / SSY,2);
                
        % R^2 for 6.25:1

        y_bar = sum(exp_data(:,2) / length(exp_data(:,2)));
        SSY = sum((exp_data(:,2)  - y_bar).^2);
        SSE = sum((exp_data(:,2)-sim_data(:,2)).^2);
        R_squared(2) = round((SSY - SSE) / SSY,2);
        
        % R^2 for 1:1

        y_bar = sum(exp_data(:,3) / length(exp_data(:,3)));
        SSY = sum((exp_data(:,3)  - y_bar).^2);
        SSE = sum((exp_data(:,3)-sim_data(:,3)).^2);
        R_squared(3) = round((SSY - SSE) / SSY,2);

    % Calculate RMSE

        % R^2 for 25:1

        SSE_tot = sum((exp_data(:,1)-sim_data(:,1)).^2) + ...
                    sum((exp_data(:,2)-sim_data(:,2)).^2) + ...
                    sum((exp_data(:,3)-sim_data(:,3)).^2);
        M_SSE = SSE_tot / (length(exp_data(:,1)) + length(exp_data(:,2)) + length(exp_data(:,3)));
        RMSE = M_SSE^0.5;

    % Create and save monte carlo parameter distribution
    
        % Save the plot in the result folder
        saveas(fig_mc, fullfile(result_folder, sprintf('monte_carlo.png')));
        close(fig_mc); % Close the figure to avoid cluttering memory           
            
    % Save calibration plot for effectors

        % Create a new figure for the current plot without displaying it
        fig_ec = plots_calibration_effector(CI_eff,cells_eff,q);
    
        % Save the plot in the result folder
        saveas(fig_ec, fullfile(result_folder, sprintf('effector_calibration.png')));
        close(fig_ec); % Close the figure to avoid cluttering memory
       
    % Save calibration plot for targets

        fig_tc = plots_calibration_target(CI_tar,cells_tar,slope,intercept);
    
        % Save the plot in the result folder
        saveas(fig_tc, fullfile(result_folder, sprintf('target_calibration.png')));
        close(fig_tc); % Close the figure to avoid cluttering memory

    % Create a CSV file from the sorted parameters and losses
    disp("Creating csv file...")
    csv_fullpath = fullfile(result_folder, csv_name);
    param_est_save_results.create_csv_from_data(sorted_pars, sorted_losses, init_pars, ...
        p_names, p_log, effector0, c0, method, seed, csv_fullpath);  
    
    % Notify user of completion
    fprintf('Plots saved to %s\n', result_folder);
    fprintf('CSV file created: %s\n', csv_fullpath);
end


%% Save parameter estimation data to file
    
function create_csv_from_data(all_pars, all_losses, param_init_guess, ...
        parameter_names, p_log, effector0, c0, sampling_method, seed, output_file)
    % Function to create a CSV from parameter estimation results.
    %
    % Inputs:
    %   all_pars         - Optimized parameter values for all iterations
    %   all_losses       - Corresponding loss values for all iterations
    %   param_init_guess - Initial guess for parameters
    %   parameter_names  - Names of the parameters
    %   p_log            - Log fit parameters (a,b)
    %   effector0        - Initial concentrations for effector data
    %   c0               - Initial concentrations for c0
    %   sampling_method  - Method used for parameter estimation ('uniform', 'gaussian', 'grid')
    %   seed             - Seed used for random number generation
    %   output_file      - Path to save the CSV file
    
    % Extract the number of iterations and parameters
    num_iterations = size(all_pars, 1);
    csv_data = cell(num_iterations, 3);  % Adjusted to have 3 columns
    
    % Preallocate for parameter names string
    param_names_str = strcat('[', strjoin(parameter_names, ' '), ']');

    % Preallocate strings for effector0 and c0 vectors
    effector0_str = strcat('[', num2str(effector0(:)'), ']');
    c0_str = strcat('[', num2str(c0(:)'), ']');

    % Print introductory lines at the beginning of the CSV
    fileID = fopen(output_file, 'w');
    fprintf(fileID, 'Model Specifications\n');
    fprintf(fileID, 'Parameter names: %s\n', param_names_str);
    fprintf(fileID, 'Log fit parameters (a b): [%g %g]\n', p_log(1), p_log(2));
    fprintf(fileID, 'Effector0: %s\n', effector0_str);
    fprintf(fileID, 'C0: %s\n', c0_str);
    fprintf(fileID, 'Method: %s\n', sampling_method);
    fprintf(fileID, 'Seed: %d\n\n', seed);
    
    % Populate the CSV data with iteration results
    for i = 1:num_iterations
        % Process initial guess for parameters
        param_init_str = sprintf('[%s]', num2str(param_init_guess(i, :)));

        % Process minimized parameters
        minimized_parameters = sprintf('[%s]', num2str(all_pars(i, :)));

        % Fill the CSV data for each iteration
        csv_data{i, 1} = num2str(all_losses(i));       % Loss
        csv_data{i, 2} = minimized_parameters;         % Minimized parameters
        csv_data{i, 3} = param_init_str;               % Initial guess
    end

    % Create table and write to CSV (including headers)
    col_names = {'loss', 'minimized_parameters', 'param_init_guess'};  % Adjusted to match data order
    csv_table = cell2table(csv_data, 'VariableNames', col_names);
    
    % Write the table and ensure it writes the column headers
    writetable(csv_table, output_file, 'WriteVariableNames', true, 'WriteMode', 'append');

    fclose(fileID);
end



    end
end