%-----------------------------------------------------------------------------
% xCellignce data fitting - target cells (raji)
% Created by: Viren Shah
% Last update date: 6/8/2023
% Prompt: target cell data fit
%-----------------------------------------------------------------------------%%

run("target_cell_data.m")

%% set variables

target_growth_rate_est = 0.05;      % Tumor growth rate (hours^-1)
carrying_capacity_est = 1e6;        % Tumor carrying capacity (cells)

%% fmincon - logfit for target growth

% options = optimset('fmincon');
% options = optimset(options,...
%     'TolFun',1e-6,...
%     'TolX',1e-8,...
%     'Display','iter',...
%     'Maxiter',500,...
%     'MaxFunEvals',5000);
% options = optimset(options,...
%     'UseParallel','always');

options = [];    % We want to use only default optimizer options

%Tumor
a = target_growth_rate_est;	            %Target growth rate (hours^-1)
inv_b_1 = carrying_capacity_est;	    %Taret carrying capacity (cells)
b_1 = 1/inv_b_1;
 
parameters0 = [a b_1]; % a b
lb = [0.005 3.3e-6]; % lower bounds
ub = [0.500 1.5e-5]; % upper bounds

time_data_fit = time_data(1:63);
mean_target_c_1_fit = mean_target_c_1(1:63);
 
[parameters_logfit,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@sys_err_log,parameters0,[],[],[],[],lb,ub,[],options,...
    time_data_fit, target_ic, mean_target_c_1_fit);

% fmincon output
% disp('Estimated parameter values'); disp(parameters_logfit);
% disp('Function at optimization termination'); disp(fval);
% disp('Exitflag at optimization termination'); disp(exitflag);
% disp('Output at optimization termination'); disp(output);
% disp('Lambda at optimization termination'); disp(lambda);
% disp('Gradient at optimization termination'); disp(grad);
% disp('Hessian at optimization termination'); disp(hessian);
                  
% model fit
[t_fit_log,target_fit_log] = ...
    ode45(@logistic_growth, time_data, target_ic, [], parameters_logfit);
target_fit_CI = (target_fit_log - b) ./ m;
residuals_logfit = (mean_target_c - target_fit_log);

