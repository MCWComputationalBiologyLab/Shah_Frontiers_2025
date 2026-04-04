classdef configuration
    methods(Static)

    function [st_index, columns_target, columns_cytolysis, columns_25_1, ...
            columns_25, columns_6_1, columns_6, columns_1_1, columns_1, ...
            target_ic, num_simulations, te, par_a, par_b, par_kappa, par_k, ...
            par_kc, par_kp2, par_kd, par_l, par_s, par_n] ...
            = extract_file_data(filename, excelData)
        % EXTRACT_FILE_DATA Extracts relevant data for a specific file from loaded Excel data.
        %
        % Parameters:
        %   filename (string) - Name of the Excel file to match in the dataset
        %   excelData (cell array) - Cell array of data loaded from the main Excel file
        %
        % Returns:
        %   st_index (int) - Start index for the file
        %   targetData (1x3 array) - Target cell growth data columns
        %   cytolysisData (1x3 array) - Target cell cytolysis data columns
        %   et_25_1_Data (1x3 array) - E:T 25:1 data columns
    
        % Extract file paths from the 'filepath' row (row 3 in your Excel file)
        filePaths = strrep(string(excelData(3, 2:end)), "'", "");

        filename = strrep(filename, '"', '');
        
        % Find the column index that matches the given filename
        colIndex = find(strcmp(filename, filePaths));
    
        if isempty(colIndex)
            error('Filename %s not found in the main Excel file.', filename);
        end  
        
        % Extract the required variables based on the matched column
        st_index = cell2mat(excelData(4, colIndex + 1)); % Row 4 for start index
        
        % Extract specific data columns for target, cytolysis, and E:T 25:1 data
        columns_target = cell2mat(excelData(5:7, colIndex + 1));  % Rows 5-7 for target data
        columns_cytolysis = cell2mat(excelData(8:10, colIndex + 1));  % Rows 8-10 for cytolysis data
        columns_25_1 = cell2mat(excelData(11:13, colIndex + 1));  % Rows 11-13 for E:T 25:1 data
        columns_25 = cell2mat(excelData(14:16, colIndex + 1));
        columns_6_1 = cell2mat(excelData(17:19, colIndex + 1));
        columns_6 = cell2mat(excelData(20:22, colIndex + 1));
        columns_1_1 = cell2mat(excelData(23:25, colIndex + 1));
        columns_1 = cell2mat(excelData(26:28, colIndex + 1));
        target_ic = cell2mat(excelData(29, colIndex + 1));
        num_simulations = cell2mat(excelData(30, colIndex + 1));
        te = cell2mat(excelData(31, colIndex + 1));
        par_a = cell2mat(excelData(32, colIndex + 1));
        par_b = cell2mat(excelData(33, colIndex + 1));
        par_kappa = cell2mat(excelData(34, colIndex + 1));
        par_k = cell2mat(excelData(35, colIndex + 1));
        par_kc = cell2mat(excelData(36, colIndex + 1));
        par_kp2 = cell2mat(excelData(37, colIndex + 1));
        par_kd = cell2mat(excelData(38, colIndex + 1));
        par_l = cell2mat(excelData(39, colIndex + 1));
        par_s = cell2mat(excelData(40, colIndex + 1));
        par_n = cell2mat(excelData(41, colIndex + 1));
    end

    function config = buildSimulationConfig(rootDir, config_file)
        % pathToConfig - path to excel file with simulation information

        % Load the main Excel file into a cell array
        excelData = readcell(fullfile(rootDir, config_file));
    
        % Initialize an empty structure to store data for each file
        config = struct();
    
        % Get a list of folders in the root directory
        folderList = dir(rootDir);
        folderList = folderList([folderList.isdir]); % Filter to only directories
    
        % Loop through each folder
        for k = 1:length(folderList)
            folder = folderList(k).name;
    
            % Skip '.' and '..' folders and folders that don't start with 'p'
            if strcmp(folder, '.') || strcmp(folder, '..') || folder(1) ~= 'p'
                continue;
            end
    
            % Construct the full path to the folder
            folderPath = fullfile(rootDir, folder);
    
            % Get a list of Excel files in this folder (assuming one Excel file per folder)
            fileList = dir(fullfile(folderPath, '*.xlsx'));
    
            % Loop through each Excel file in the folder
            for j = 1:length(fileList)
                dataFileName = fileList(j).name;

                % Skip any files that don't start with RAJI
                if ~strncmp(dataFileName, 'RAJI',4)
                    continue;
                end

                % Clean name to create data struct path
                fName = matlab.lang.makeValidName(dataFileName);
    
                % Extract data for the specific file
                [st_index, columns_target, columns_cytolysis, columns_25_1, ...
                columns_25, columns_6_1, columns_6, columns_1_1, columns_1, ...
                target_ic, num_simulations, te, par_a, par_b, par_kappa, par_k, ...
                par_kc, par_kp2, par_kd, par_l, par_s, par_n] ...
                = configuration.extract_file_data(dataFileName, excelData);
                
                % Load data file
                %data = readmatrix(fullfile(folderPath, dataFileName));
                
                data = configuration.datasetimport(folderPath, dataFileName);

                % Store extracted data in the structure
                config.(folder).(fName).st_index = st_index;
                config.(folder).(fName).columns_target = columns_target';
                config.(folder).(fName).columns_cytolysis = columns_cytolysis';
                config.(folder).(fName).columns_25_1 = columns_25_1';
                config.(folder).(fName).columns_25 = columns_25';
                config.(folder).(fName).columns_6_1 = columns_6_1';
                config.(folder).(fName).columns_6 = columns_6';
                config.(folder).(fName).columns_1_1 = columns_1_1';
                config.(folder).(fName).columns_1 = columns_1';
                config.(folder).(fName).target_ic = target_ic';
                config.(folder).(fName).num_simulations = num_simulations;
                config.(folder).(fName).te = te;
                config.(folder).(fName).par_a = par_a;
                config.(folder).(fName).par_b = par_b;
                config.(folder).(fName).par_kappa = par_kappa;
                config.(folder).(fName).par_k = par_k;
                config.(folder).(fName).par_kc = par_kc;
                config.(folder).(fName).par_kp2 = par_kp2;
                config.(folder).(fName).par_kd = par_kd;
                config.(folder).(fName).par_l = par_l;
                config.(folder).(fName).par_s = par_s;
                config.(folder).(fName).par_n = par_n;
                config.(folder).(fName).data = data;

            end
        end
    end


    function data = datasetimport(folderPath, dataFileName)

        %% Import datafile

        dataset = readtable(fullfile(folderPath, dataFileName), 'Sheet', 'Well Graph');
        data = table2array(dataset);
        disp(fullfile(folderPath, dataFileName))

        % Check for NaNs in the first 5 columns
        nanRows = any(isnan(data(:, 1:5)), 2);

        % Remove rows with NaNs in the first 5 columns
        data_clean = data(~nanRows, :);

        % Remove duplicate rows
        data_clean = unique(data_clean, 'rows');

        %% Set data

        data = data_clean;

    end

    end
end
