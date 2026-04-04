%-----------------------------------------------------------------------------
% xCellignce data processing - effector cells (CART)
% Created by: Viren Shah
% Last update date: 6/8/2023
% Prompt: effector cell dataset (run after targets)
%-----------------------------------------------------------------------------%%

%% Effector CI to cell conversion - decay

Ein_decay = Ein+1;
CI_25_decay = mean_e25(Ein_decay);
CI_6_decay = mean_e6(Ein_decay);
CI_1_decay = mean_e1(Ein_decay);
CI_0_decay = mean_e1(Ein_decay-4);

CI_eff_decay = [CI_0_decay; CI_1_decay; CI_6_decay; CI_25_decay;];
cells_eff = [0; ...
         (target_ic); ...
         (target_ic*6.25); ...
         (target_ic*25)];

qfit_decay = polyfit(cells_eff, CI_eff_decay, 1);
cells_range_decay = (0:1000:2e6);
CI_calc_decay = polyval(qfit_decay,cells_range_decay);
CI_calc_decay_1 = polyval(qfit_decay,cells_eff);

        y_bar = sum(CI_eff_decay / length(CI_eff_decay));
        SSY = sum((CI_eff_decay  - y_bar).^2);
        SSE = sum((CI_eff_decay-CI_calc_decay_1).^2);
        R_squared_eff = round((SSY - SSE) / SSY,2);

%% fmincon - effector decay assays

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
kd = 1e-3;	            %Tumor growth rate (hours^-1)

parameters0 = [kd]; % kd
lb = [0.000001]; % lower bounds
ub = [0.050000]; % upper bounds

effector_ic = [cells_eff(2) 
               cells_eff(3) 
               cells_eff(4)];
mean_effectors_efit = [mean_e1(Ein_decay:end),...
                       mean_e6(Ein_decay:end),...
                       mean_e25(Ein_decay:end)];
std_effectors_efit = [std_e1(Ein_decay:end),...
                       std_e6(Ein_decay:end),...
                       std_e25(Ein_decay:end)];
Tdata_efit = time_data(Ein_decay:end); 
[e_kd,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@sys_err_decay,parameters0,[],[],[],[],lb,ub,[],options,...
    Tdata_efit, effector_ic, mean_effectors_efit, qfit_decay);

% fmincon output
% disp('Estimated parameter values'); disp(parameters_logfit);
% disp('Function at optimization termination'); disp(fval);
% disp('Exitflag at optimization termination'); disp(exitflag);
% disp('Output at optimization termination'); disp(output);
% disp('Lambda at optimization termination'); disp(lambda);
% disp('Gradient at optimization termination'); disp(grad);
% disp('Hessian at optimization termination'); disp(hessian);
                  
    % model fit
    [t_fit_decay,effector_fit_decay_1] = ...
        ode45(@effector_decay, Tdata_efit, effector_ic(1), [], e_kd);
    effector_fit_CI_1 = polyval(qfit_decay,effector_fit_decay_1);
    [t_fit_decay,effector_fit_decay_6] = ...
        ode45(@effector_decay, Tdata_efit, effector_ic(2), [], e_kd);
    effector_fit_CI_6 = polyval(qfit_decay,effector_fit_decay_6);
    [t_fit_decay,effector_fit_decay_25] = ...
        ode45(@effector_decay, Tdata_efit, effector_ic(3), [], e_kd);
    effector_fit_CI_25 = polyval(qfit_decay,effector_fit_decay_25);