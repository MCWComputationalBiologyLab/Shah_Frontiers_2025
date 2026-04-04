clc; close all; clearvars;

%% Set Model Run Parameters

% Use parallel processing
use_parallel = true;

% sampling_method: uniform, gaussian, grid
% "uniform" - randomly pull points between lb and ub uniformly
% "gaussian" - randomly pull points using a parameter vector as the mean
    % value and the bounds to establish a standard deviation (bell curve)
% "grid" - exhaustively try every combination of parameters in a grid.  This
    % method is very intensive the more grid points and parameters in model.
    % (ignores num_simulations)

sampling_method = "gaussian";

% Only used if method "grid" is used

num_grid_points = 5;

% Set rng value for reproducibility

seed = 1;
rng(seed);

% Parameters for selecting results

top_n = 5; % Number of top loss minimizations to plot and save 
columns = [1,2,3]; % Columns from result/exp matrices to compare

%% Set directory

% Define the root directory 

rootDir = pwd;

%% Create output CSVs

% Create result CSV, delete previous if exists

resultCSV = fullfile(rootDir, 'results_summary_p99_uniform_te_m1.csv');

if exist(resultCSV, 'file') == 2 
    delete(resultCSV);  % Remove old summary if exists
end
patCSV = fullfile(rootDir,'patCSV.csv');

if exist(patCSV, 'file') == 2
    delete(patCSV);  % Remove patient file if exists
end
timeCSV = fullfile(rootDir,'model_time.csv');

if exist(timeCSV, 'file') == 2
    delete(timeCSV);  % Remove time if exists
end
tumorCSV = fullfile(rootDir,'tumor.csv');

if exist(tumorCSV, 'file') == 2
    delete(tumorCSV);  % Remove tumor only data if exists
end
e25t1CSV = fullfile(rootDir,'e25t1.csv');

if exist(e25t1CSV, 'file') == 2
    delete(e25t1CSV);  % Remove effector data if exists
end
e6t1CSV = fullfile(rootDir,'e6t1.csv');

if exist(e6t1CSV, 'file') == 2
    delete(e6t1CSV);  % Remove effector data if exists
end
e1t1CSV = fullfile(rootDir,'e1t1.csv');

if exist(e1t1CSV, 'file') == 2
    delete(e1t1CSV);  % Remove effector data if exists
end
cytolysisCSV = fullfile(rootDir,'cytolysis.csv');

if exist(cytolysisCSV, 'file') == 2
    delete(cytolysisCSV);  % Remove cytolysis only data if exists
end
e25CSV = fullfile(rootDir,'e25.csv');

if exist(e25CSV, 'file') == 2
    delete(e25CSV);  % Remove effector data if exists
end
e6CSV = fullfile(rootDir,'e6.csv');

if exist(e6CSV, 'file') == 2
    delete(e6CSV);  % Remove effector data if exists
end
e1CSV = fullfile(rootDir,'e1.csv');

if exist(e1CSV, 'file') == 2
    delete(e1CSV);  % Remove effector data if exists
end
e25t1_sim_CI_CSV = fullfile(rootDir,'e25t1_sim_CI.csv');

if exist(e25t1_sim_CI_CSV, 'file') == 2
    delete(e25t1_sim_CI_CSV);  % Remove sim data if exists
end
e6t1_sim_CI_CSV = fullfile(rootDir,'e6t1_sim_CI.csv');

if exist(e6t1_sim_CI_CSV, 'file') == 2
    delete(e6t1_sim_CI_CSV);  % Remove sim data if exists
end
e1t1_sim_CI_CSV = fullfile(rootDir,'e1t1_sim_CI.csv');

if exist(e1t1_sim_CI_CSV, 'file') == 2
    delete(e1t1_sim_CI_CSV);  % Remove sim data if exists
end
target_sim_CI_CSV = fullfile(rootDir,'target_sim_CI.csv');

if exist(target_sim_CI_CSV, 'file') == 2
    delete(target_sim_CI_CSV);  % Remove sim data if exists
end
Ein_CSV = fullfile(rootDir,'Ein.csv');

if exist(Ein_CSV, 'file') == 2
    delete(Ein_CSV);  % Remove Ein data if exists
end

% Initialize the result CSV file
file1 = fopen(resultCSV, 'w');
fprintf(file1, 'Folder, loss, a, b, Kappa, k, kc, kp2, kd, l, s, R2_25, R2_6, R2_1, RMSE R2_EC\n'); 
fclose(file1);

%% Set configuration file 

config_file = "simulation_input.xlsx";
new_config_file = 1;

if new_config_file == 0
    if exist("config_struct.mat", 'file')
    load("config_struct.mat");
    else
    config_struct = configuration.buildSimulationConfig(rootDir, config_file);
    save("config_struct.mat");
    end
else
    config_struct = configuration.buildSimulationConfig(rootDir, config_file);
    save("config_struct.mat");
end

%% Set file variables

for folder = fieldnames(config_struct)'

    folderName = folder{1};
    disp(strcat('Processing :~', folderName));

    % Should only be one here (TODO Concat)
    for file = fieldnames(config_struct.(folderName))'
    disp(strcat('Processing :~', file));
    fileName = file{1};
        
        % Extract the variables from the current structure entry

        st_index = config_struct.(folderName).(fileName).st_index;
        columns_target = config_struct.(folderName).(fileName).columns_target;
        columns_cytolysis = config_struct.(folderName).(fileName).columns_cytolysis;
        columns_25_1 = config_struct.(folderName).(fileName).columns_25_1;
        columns_25 = config_struct.(folderName).(fileName).columns_25;
        columns_6_1 = config_struct.(folderName).(fileName).columns_6_1;
        columns_6 = config_struct.(folderName).(fileName).columns_6;
        columns_1_1 = config_struct.(folderName).(fileName).columns_1_1;
        columns_1 = config_struct.(folderName).(fileName).columns_1;
        target_ic = config_struct.(folderName).(fileName).target_ic;
        num_simulations = config_struct.(folderName).(fileName).num_simulations;
        te = config_struct.(folderName).(fileName).te;
        par_a = config_struct.(folderName).(fileName).par_a;
        par_b = config_struct.(folderName).(fileName).par_b;
        par_kappa = config_struct.(folderName).(fileName).par_kappa;
        par_k = config_struct.(folderName).(fileName).par_k;
        par_kc = config_struct.(folderName).(fileName).par_kc;
        par_kp2 = config_struct.(folderName).(fileName).par_kp2;
        par_kd = config_struct.(folderName).(fileName).par_kd;
        par_l = config_struct.(folderName).(fileName).par_l;
        par_s = config_struct.(folderName).(fileName).par_s;
        par_n = config_struct.(folderName).(fileName).par_n;
        
        %number of simulations override
        num_simulations = 1000;
    
        data = config_struct.(folderName).(fileName).data;
        
        % clear variables for next run

        clearvars -except folderName st_index columns_target ...
                          columns_cytolysis columns_25_1 columns_25 ...
                          columns_6_1 columns_6 columns_1_1 columns_1 ...
                          target_ic num_simulations te par_a par_b par_kappa ...
                          par_k par_kc par_kp2 par_kd par_l par_s par_n ...
                          data rootDir config_file ...
                          resultCSV timeCSV tumorCSV ...
                          e25t1CSV e6t1CSV e1t1CSV patCSV ...
                          cytolysisCSV e25CSV e6CSV e1CSV ...
                          e25t1_sim_CI_CSV e6t1_sim_CI_CSV e1t1_sim_CI_CSV ...
                          target_sim_CI_CSV Ein_CSV ...
                          config_struct use_parallel sampling_method ...
                          num_grid_points seed columns top_n R_squared_eff

        % Create a folder-specific results CSV file path

        save_dir = fullfile(rootDir, folderName);
        csv_name = "results.csv";

        % Run dataload and parameterization for targets (tumor) only

        run("target_cell_driver.m")
    
        % Run dataload and parameterization for effector-target (CAR T-cell - tumor) cytotoxicity
        % experiments

        run("effector_cell_driver.m")

        % Save results

        [target_fit_CI,...
            sim_data, sim_data_tc, sim_data_ec,...
            sim_data_utc, sim_data_fc, sim_data_fp, R_squared, RMSE] ...
         = param_est_save_results.saveResults(all_losses, all_pars, init_pars, tspan, c0, p_log, ...
            effector0, tc0, mean_target_CI, std_target_CI, mean_target_c, target_fit_CI, target_fit_log, ...
            mean_targets, std_targets, columns, cells_tar, CI_tar, cells_eff, CI_eff,...
            slope, intercept, qfit, Ein, et, top_n, ...
            save_dir, p_names, csv_name, seed, sampling_method, false);

        % From this above we get all_pars, all_losses from the parest
        % Sort the all_losses array in ascending order

        [sorted_losses, sort_idx] = sort(all_losses);
        
        % Reorder all_pars to match the sorted all_losses

        sorted_pars = all_pars(sort_idx, :);

        % Get the top value in sorted_losses and the corresponding vector from sorted_pars

        loss = sorted_losses(1);  % Top value in sorted_losses
        top_par = sorted_pars(1, :);  % Assuming sorted_pars is a matrix where each row is a parameter vector
        
        % Open the CSV file in append mode to add the results     
        % Write the results to the CSV file
        % Close the file

        % Output SSR parameters

        file1 = fopen(resultCSV, 'a');
        fprintf(file1, '%s, %.4f', folderName, loss);  % Folder name and loss value
        fprintf(file1, ', %.4e', parameters_logfit);
        fprintf(file1, ', %.4f', top_par);  % Write each parameter in the vector
        fprintf(file1, ', %.4f', R_squared);  % Write each parameter in the vector
        fprintf(file1, ', %.4f', RMSE);  % Write each parameter in the vector
        fprintf(file1, ', %.4f', R_squared_eff);  % Write each parameter in the vector
        fprintf(file1, '\n');  % New line after each row
        fclose(file1);

        % Output time
        file2 = fopen(timeCSV, 'a');
        fprintf(file2, '%.4f,', tspan);
        fprintf(file2, '\n');  % New line after each row
        fclose(file2);

        % Output tumor only data      
        file3 = fopen(tumorCSV, 'a');
        fprintf(file3, '%.4f,', mean_target_CI);
        fprintf(file3, '\n');  % New line after each row
        fclose(file3);
        
        % Output effector:tumor (25:1) data
        file4 = fopen(e25t1CSV, 'a');
        fprintf(file4, '%.4f,', mean_e25t1);
        fprintf(file4, '\n');  % New line after each row
        fclose(file4);
        
        % Output effector:tumor (6:1) data    
        file5 = fopen(e6t1CSV, 'a');
        fprintf(file5, '%.4f,', mean_e6t1);
        fprintf(file5, '\n');  % New line after each row
        fclose(file5);
        
        % Output effector:tumor (1:1) data    
        file6 = fopen(e1t1CSV, 'a');         
        fprintf(file6, '%.4f,', mean_e1t1);
        fprintf(file6, '\n');  % New line after each row
        fclose(file6);

        % Output patient
        file7 = fopen(patCSV, 'a');
        fprintf(file7, '%s', folderName);
        fprintf(file7, '\n');  % New line after each row
        fclose(file7);

        % Output cytolysis
        file8 = fopen(cytolysisCSV, 'a');
        fprintf(file8, '%.4f,', mean_cytolysis_CI);
        fprintf(file8, '\n');  % New line after each row
        fclose(file8);

        % Output e25
        file9 = fopen(e25CSV, 'a');
        fprintf(file9, '%.4f,', mean_e25);
        fprintf(file9, '\n');  % New line after each row
        fclose(file9);

        % Output e6
        file10 = fopen(e6CSV, 'a');
        fprintf(file10, '%.4f,', mean_e6);
        fprintf(file10, '\n');  % New line after each row
        fclose(file10);

        % Output e1
        file11 = fopen(e1CSV, 'a');
        fprintf(file11, '%.4f,', mean_e1);
        fprintf(file11, '\n');  % New line after each row
        fclose(file11);

        % Output e25t1 sim
        file12 = fopen(e25t1_sim_CI_CSV, 'a');
        fprintf(file12, '%.4f,', sim_data(:,1));
        fprintf(file12, '\n');  % New line after each row
        fclose(file12);

        % Output e6t1 sim
        file13 = fopen(e6t1_sim_CI_CSV, 'a');
        fprintf(file13, '%.4f,', sim_data(:,2));
        fprintf(file13, '\n');  % New line after each row
        fclose(file13);

        % Output e1t1 sim
        file14 = fopen(e1t1_sim_CI_CSV, 'a');
        fprintf(file14, '%.4f,', sim_data(:,3));
        fprintf(file14, '\n');  % New line after each row
        fclose(file14);

        % Output e1t1 sim
        file15 = fopen(target_sim_CI_CSV, 'a');
        fprintf(file15, '%.4f,', target_fit_CI);
        fprintf(file15, '\n');  % New line after each row
        fclose(file15);

        % Output Ein
        file16 = fopen(Ein_CSV, 'a');
        fprintf(file16, '%.4f,', Ein);
        fprintf(file16, '\n');  % New line after each row
        fclose(file16);
        
    end
end
