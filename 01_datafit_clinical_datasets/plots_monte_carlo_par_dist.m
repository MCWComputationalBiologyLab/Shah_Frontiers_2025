function f = plots_monte_carlo_par_dist(sorted_pars,sorted_losses)

f = figure("Visible","off");
set(gcf,'Units','inches','Position',[0.5,0.5,16,4])

parameters_order(:,1) = sorted_pars(:,3);
parameters_order(:,2) = sorted_pars(:,7);
parameters_order(:,3) = sorted_pars(:,6);
parameters_order(:,4) = sorted_pars(:,4);
parameters_order(:,5) = sorted_pars(:,2);
parameters_order(:,6) = sorted_pars(:,1);
parameters_order(:,7) = sorted_losses;
sorted_pars = parameters_order;
% ylim_kc = [-0.05 0.37]; ylim_s = [-0.90 4]; ylim_l = [-0.90 4]; ylim_kp2 = [-0.05 0.37];
% ylim_k = [1e2 4e5]; ylim_kappa = [1e4 5e6];
% xlim_par = [ylim_kc; ylim_s; ylim_l; ylim_kp2; ylim_k; ylim_kappa];
par_names = {'k_c','K_{mr}','n_{ }','k_{p2}','K_{mp}','C_E','SSR'};
x_exp_par = [0 0 0 0 3 5 0];

t = tiledlayout(1,7);
for i = 1:width(sorted_pars)
nexttile
param = categorical(cellstr(par_names(i)));
histogram(sorted_pars(:,i));

    ax = gca;
    ax.XAxis.Exponent = x_exp_par(i);
%    xlim(xlim_par(i,:));

xlabel(param)
set(gca,'Fontsize',18); 

end

end