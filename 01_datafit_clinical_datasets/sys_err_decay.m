function SSR = sys_err_decay(parameters, tspan, effector_ic, mean_effectors_efit, q)

[t,effector_fit1] = ode45(@effector_decay, tspan, effector_ic(1), [], parameters);
[t,effector_fit6] = ode45(@effector_decay, tspan, effector_ic(2), [], parameters);
[t,effector_fit25] = ode45(@effector_decay, tspan, effector_ic(3), [], parameters);

effector_fit1_CI = polyval(q,effector_fit1);
effector_fit6_CI = polyval(q,effector_fit6);
effector_fit25_CI = polyval(q,effector_fit25);

SSR_pts = ((mean_effectors_efit(:,1) - effector_fit1_CI).^2) + ...
          ((mean_effectors_efit(:,2) - effector_fit6_CI).^2) + ...
          ((mean_effectors_efit(:,3) - effector_fit25_CI).^2);

SSR = sum(SSR_pts,'all');

end