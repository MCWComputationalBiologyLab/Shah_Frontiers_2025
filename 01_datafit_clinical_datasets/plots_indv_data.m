function f = plots_indv_data(time, mean_target_CI, std_target_CI, mean_target_c, target_fit_CI, target_fit_log,  ...
                    exp_data, exp_data_std, sim_data, sim_data_tc, sim_data_ec, sim_data_utc, sim_data_fc, sim_data_fp, columns, Ein)
   
    % Set figure variables
    f = figure("Visible","off");

    set(gcf,'Units','inches','Position',[0.5,0.5,15,15]);
    
    s = 0;
    m = 4;
    n = 3;
    
    % Targets only data fit
    s = s+1;
    subplot(m,n,s)
    plot(time,mean_target_CI,':r','Linewidth',2.5); hold on;
    patch_time = [time' flip(time')];
    patch_array = [(mean_target_CI + std_target_CI)' flip((mean_target_CI - std_target_CI)')];
    patch(patch_time, patch_array, [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    plot(time,target_fit_CI,'-r','Linewidth',2.5);
    title('Target Model Fit')
    ylim([-0.5 4]);
    yticks([0 1 2 3 4]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('CI');
    set(gca,'Fontsize',20); 
    H = legend('Data', 'SD', 'Model','Location','south');
    set(H,'FontSize',16)
    legend boxoff

    % Target only conversion into cells
    s = s+1;
    subplot(m,n,s)
    plot(time,mean_target_c,':r','Linewidth',2.5); hold on;
    plot(time,target_fit_log,'-r','Linewidth',2.5);
    title('Target (Cells)')
    ylim([0 10000000]);
    yticks([10^3 10^4 10^5 10^6]);
    ylim([10^3 10^6]);
    xlim([0 80]);
    set(gca, 'YScale', 'log');
    set(gca,'Fontsize',20); 
    xlabel('Time (hours)'); ylabel('# Cells');
    H = legend('Converted Cells', 'Model','Location','south');
    set(H,'FontSize',16)
    legend boxoff

    % CI data fit
    s = s+2;
    subplot(m,n,s)
    plot(time,exp_data(:,1),':b','Linewidth',2.5); hold on;
    patch_time = [time' flip(time')];
    patch_array = [(exp_data(:,1) + exp_data_std(:,1))' flip((exp_data(:,1) - exp_data_std(:,1))')];
    patch(patch_time, patch_array, [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    plot(time,sim_data(:,1),'-b','Linewidth',2.5); hold on;
    title('E:T (25:1) Model Fit')
    ylim([-0.5 4]);
    yticks([0 1 2 3 4]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('CI');
    set(gca,'Fontsize',20);     
    H = legend('Data', 'SD', 'Model','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff
    
    s = s+1;
    subplot(m,n,s)
    plot(time,exp_data(:,2),':m','Linewidth',2.5); hold on;
    patch_time = [time' flip(time')];
    patch_array = [(exp_data(:,2) + exp_data_std(:,1))' flip((exp_data(:,2) - exp_data_std(:,2))')];
    patch(patch_time, patch_array, [1 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    plot(time,sim_data(:,2),'-m','Linewidth',2.5); hold on;
    title('E:T (6.25:1) Model Fit')
    ylim([-0.5 4]);
    yticks([0 1 2 3 4]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('CI');
    set(gca,'Fontsize',20);     
    H = legend('Data', 'SD', 'Model','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff
    
    s = s+1;
    subplot(m,n,s)
    plot(time,exp_data(:,3),':g','Linewidth',2.5); hold on;
    patch_time = [time' flip(time')];
    patch_array = [(exp_data(:,3) + exp_data_std(:,3))' flip((exp_data(:,3) - exp_data_std(:,3))')];
    patch(patch_time, patch_array, [0 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    plot(time,sim_data(:,3),'-g','Linewidth',2.5); hold on;
    title('E:T (1:1) Model Fit')
    ylim([-0.5 4]);
    yticks([0 1 2 3 4]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('CI');
    set(gca,'Fontsize',20);     
    H = legend('Data', 'SD', 'Model','Location','southeast');
    set(H,'FontSize',16)
    legend boxoff
    
    % Results in cells
    s = s+1;
    subplot(m,n,s)
    plot(time,sim_data_tc(:,1),'-b','Linewidth',2.5); hold on;
    plot(time,sim_data_ec(:,1),'--b','Linewidth',2.5); hold on;
    plot(time,sim_data_utc(:,1),':b','Linewidth',2.5); 
    set(gca, 'YScale', 'log')
    ylim([10^3 10^7]);
    yticks([10^3 10^4 10^5 10^6 10^7]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('# Cells');
    title('E:T (25:1) Cells')
    set(gca,'Fontsize',20); 
    H = legend('Target','CART','UTC','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff
    
    s = s+1;
    subplot(m,n,s)
    plot(time,sim_data_tc(:,2),'-m','Linewidth',2.5); hold on;
    plot(time,sim_data_ec(:,2),'--m','Linewidth',2.5); hold on;
    plot(time,sim_data_utc(:,2),':m','Linewidth',2.5); 
    set(gca, 'YScale', 'log')
    ylim([10^3 10^7]);
    yticks([10^3 10^4 10^5 10^6 10^7]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('# Cells');
    title('E:T (6.25:1) Cells')
    set(gca,'Fontsize',20); 
    H = legend('Target','CART','UTC','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff
    
    s = s+1;
    subplot(m,n,s)
    plot(time,sim_data_tc(:,3),'-g','Linewidth',2.5); hold on;
    plot(time,sim_data_ec(:,3),'--g','Linewidth',2.5); hold on;
    plot(time,sim_data_utc(:,3),':g','Linewidth',2.5);
    set(gca, 'YScale', 'log')
    ylim([10^3 10^7]);
    yticks([10^3 10^4 10^5 10^6 10^7]);
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('# Cells');
    title('E:T (1:1) Cells')
    set(gca,'Fontsize',20); 
    H = legend('Target','CART','UTC','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff

    % cytolysis
    s = s+1;
    subplot(m,n,s)
    plot(time(Ein:end),((1-(sim_data_tc(Ein:end,1) / sim_data_tc(Ein-1,1)))*100),'--r','Linewidth',2.5); hold on;
    plot(time(Ein:end),((1-(sim_data_tc(Ein:end,2) / sim_data_tc(Ein-1,2)))*100),':r','Linewidth',2.5); hold on;
    plot(time(Ein:end),((1-(sim_data_tc(Ein:end,3) / sim_data_tc(Ein-1,3)))*100),'.r','Linewidth',2.5); hold on;
    ylim([-25 125]);
    yticks([0 20 40 60 80 100]);
    set(gca, 'YScale', 'linear')
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('% cytolysis');
    title('Target Cytolysis (Cells)')
    set(gca,'Fontsize',20); 
    H = legend('25:1','6.25:1','1:1','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff

    % fc
    s = s+1;
    subplot(m,n,s)
    plot(time(Ein:end),sim_data_fc(Ein:end,1),'--r','Linewidth',2.5); hold on;
    plot(time(Ein:end),sim_data_fc(Ein:end,2),':r','Linewidth',2.5); hold on;
    plot(time(Ein:end),sim_data_fc(Ein:end,3),'.r','Linewidth',2.5); hold on;
    %ylim([-25 125]);
    %yticks([0 20 40 60 80 100]);
    set(gca, 'YScale', 'linear')
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('# Cells');
    title('FC')
    set(gca,'Fontsize',20); 
    H = legend('25:1','6.25:1','1:1','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff

    % fp
    s = s+1;
    subplot(m,n,s)
    plot(time(Ein:end),sim_data_fp(Ein:end,1),'--r','Linewidth',2.5); hold on;
    plot(time(Ein:end),sim_data_fp(Ein:end,2),':r','Linewidth',2.5); hold on;
    plot(time(Ein:end),sim_data_fp(Ein:end,3),'.r','Linewidth',2.5); hold on;
    %ylim([-25 125]);
    %yticks([0 20 40 60 80 100]);
    set(gca, 'YScale', 'linear')
    xlim([0 80]);
    xticks([0 20 40 60 80]);
    xlabel('Time (hours)'); ylabel('# Cells');
    title('FP')
    set(gca,'Fontsize',20); 
    H = legend('25:1','6.25:1','1:1','Location','northeast');
    set(H,'FontSize',16)
    legend boxoff
    
end