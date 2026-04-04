function SSR = sys_err_log(parameters, tspan, targets_ic, mean_targets)

[t,target_fit] = ode45(@logistic_growth, tspan, targets_ic, [], parameters);

SSR_pts = ((mean_targets - target_fit).^2);
SSR = sum(SSR_pts,'all');

end