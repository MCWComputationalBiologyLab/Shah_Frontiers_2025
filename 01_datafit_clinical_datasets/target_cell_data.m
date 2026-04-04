%-----------------------------------------------------------------------------
% xCellignce data processing - target cells (raji)
% Created by: Viren Shah
% Last update date: 6/8/2023
% Prompt: target cell dataset
%-----------------------------------------------------------------------------%%

%% identify time and calculate

time_data = data(st_index:end,1) ./ 3600;       % hours
et = size(time_data,1);                         % end time

%% identify targets, average, and save

target_CI = data(st_index:end,columns_target);              % target cell reference in CI
[std_target_CI,mean_target_CI] = std(target_CI,0,2);        % mean and std

%% target Cytolysis

cytolysis_CI = data(st_index:end,columns_cytolysis);            % target cell cytolysis reference in CI
[std_cytolysis_CI,mean_cytolysis_CI] = std(cytolysis_CI,0,2);   % mean and std

%% E25:T1

e25t1 = data(st_index:end,columns_25_1);             % A1 - A3 CI
[std_e25t1,mean_e25t1] = std(e25t1,0,2);             % mean and std

e25 = data(st_index:end,columns_25);               % A4 - A6 CI
[std_e25,mean_e25] = std(e25,0,2);          % mean and std

%% E6.25:T1

e6t1 = data(st_index:end,columns_6_1);            % B1 - B3 CI
[std_e6t1,mean_e6t1] = std(e6t1,0,2);       % mean and std

e6 = data(st_index:end,columns_6);              % B4 - B6 CI
[std_e6,mean_e6] = std(e6,0,2);             % mean and std

%% E1:T1

e1t1 = data(st_index:end,columns_1_1);            % C1 - C3 CI
[std_e1t1,mean_e1t1] = std(e1t1,0,2);       % mean and std

e1 = data(st_index:end,columns_1);              % C4 - C6 CI
[std_e1,mean_e1] = std(e1,0,2);             % mean and std

%% effector input location - algorithm / cross check

    % zero effectors
    [~, Emin] = min(mean(mean_e1(:,1) + mean_e6(:,1) + mean_e25(:,1))/3);

    %effector input
    mean_e25_m1(1) = 0;
    for i = 2:size(mean_e25,1)
    mean_e25_m1(i) = mean_e25(i) - mean_e25(i-1);
    end
    mean_e25_m1 = mean_e25_m1';
    mean_e25_de = mean_e25 + mean_e25_m1;
    [~,Ein] = max(mean_e25_de(1:30,1));

%% Dataset preprocessing

min_CI = min([mean_target_CI mean_cytolysis_CI ...
               mean_e25t1 mean_e25 ...
               mean_e6t1 mean_e6 ...
               mean_e1t1 mean_e1],[],'all');

mean_target_CI = mean_target_CI - min_CI;
mean_cytolysis_CI = mean_cytolysis_CI - min_CI;
mean_e25t1 = mean_e25t1 - min_CI;
mean_e25 = mean_e25 - min_CI;
mean_e6t1 = mean_e6t1 - min_CI;
mean_e6 = mean_e6 - min_CI;
mean_e1t1 = mean_e1t1 - min_CI;
mean_e1 = mean_e1 - min_CI;

drift = max(mean_e1(Ein) - mean_e1(Ein:end),0);
mean_e25(Ein:end) = mean_e25(Ein:end) + drift;
mean_e6(Ein:end) = mean_e6(Ein:end) + drift;
mean_e1(Ein:end) = mean_e1(Ein:end) + drift;

s_CI = 1.00;
diff_CI = s_CI - mean_target_CI(1);
mean_target_CI = mean_target_CI + diff_CI;
diff_CI = s_CI - mean_cytolysis_CI(1);
mean_cytolysis_CI = mean_cytolysis_CI + diff_CI;
diff_CI = s_CI - mean_e25t1(1);
mean_e25t1 = mean_e25t1 + diff_CI;
diff_CI = s_CI - mean_e6t1(1);
mean_e6t1 = mean_e6t1 + diff_CI;
diff_CI = s_CI - mean_e1t1(1);
mean_e1t1 = mean_e1t1 + diff_CI;

CI_min = min([mean_target_CI mean_cytolysis_CI ...
               mean_e25t1 mean_e25 ...
               mean_e6t1 mean_e6 ...
               mean_e1t1 mean_e1],[],'all');

CI_max = max([mean_target_CI mean_cytolysis_CI ...
               mean_e25t1 mean_e25 ...
               mean_e6t1 mean_e6 ...
               mean_e1t1 mean_e1],[],'all');
ls = 0.0;
us = 4.0;
mean_target_CI = rescale(mean_target_CI,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_cytolysis_CI = rescale(mean_cytolysis_CI,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_e25t1 = rescale(mean_e25t1,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_e25 = rescale(mean_e25,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_e6t1 = rescale(mean_e6t1,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_e6 = rescale(mean_e6,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_e1t1 = rescale(mean_e1t1,ls,us,"InputMin",CI_min,"InputMax",CI_max);
mean_e1 = rescale(mean_e1,ls,us,"InputMin",CI_min,"InputMax",CI_max);

s_CI = 1.00;
diff_CI = s_CI - mean_target_CI(1);
mean_target_CI = mean_target_CI + diff_CI;
diff_CI = s_CI - mean_cytolysis_CI(1);
mean_cytolysis_CI = mean_cytolysis_CI + diff_CI;
diff_CI = s_CI - mean_e25t1(1);
mean_e25t1 = mean_e25t1 + diff_CI;
diff_CI = s_CI - mean_e6t1(1);
mean_e6t1 = mean_e6t1 + diff_CI;
diff_CI = s_CI - mean_e1t1(1);
mean_e1t1 = mean_e1t1 + diff_CI;

target_growth_CI = [mean_target_CI(1:Ein-1) ...
                 mean_e25t1(1:Ein-1) ...
                 mean_e6t1(1:Ein-1) ...
                 mean_e1t1(1:Ein-1)];
[std_target_growth_CI,mean_target_growth_CI] = std(target_growth_CI,0,2);
differential = mean_target_growth_CI(end) - mean_target_CI(Ein-1);
mean_target_CI_2 = [mean_target_growth_CI; mean_target_CI(Ein:et) + differential];
mean_target_CI = mean_target_CI_2;

%% Target CI to cell conversion

CI_tar = [min(mean_cytolysis_CI); mean_target_CI(1)];
cells_tar = [0; 40000];
m = ((0 - target_ic) / (min(mean_cytolysis_CI) - mean_target_CI(1)));
b = (-1)*m*min(mean_cytolysis_CI);
mean_target_c = mean_target_CI*m + b;
mean_target_c_1(1:et) = mean_target_c(1:et);
mean_target_c_1 = mean_target_c_1';
mean_cytolysis_c = mean_cytolysis_CI*m + b;

%% Effector CI to cell conversion

CI_0 = 0;
CI_25 = mean_e25(Ein+2)-mean_e25(Ein-2);
CI_6 = mean_e6(Ein+2)-mean_e6(Ein-2);
CI_1 = mean_e1(Ein+2)-mean_e1(Ein-2);
if CI_25 > CI_6
    CI_eff = [CI_0; CI_1; CI_6; CI_25;];
else
    CI_0 = 0;
    CI_25 = mean_e25t1(Ein)-mean_e25t1(Ein-1);
    CI_6 = mean_e6t1(Ein)-mean_e6t1(Ein-1);
    CI_1 = mean_e1t1(Ein)-mean_e1t1(Ein-1);
    CI_eff = [CI_0; CI_1; CI_6; CI_25;];
end
cells_eff = [0; ...
         (target_ic); ...
         (target_ic*6.25); ...
         (target_ic*25)];
pfit = polyfit(CI_eff, cells_eff, 1);
CI_range = (0:0.01:1);
cells_calc = polyval(pfit,CI_range);
qfit = polyfit(cells_eff, CI_eff, 1);
cells_range = (0:1000:2000000);
CI_calc = polyval(qfit,cells_range);
