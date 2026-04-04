clc; close all; clearvars;

%% Load results data file

s_results = readcell("results_summary_uniform_m1.csv");

%% Load patient data file

p_data = readcell("CART_collated_data.xlsx");
p_data = p_data';

%% Load xcelligence datasets

time = readmatrix("model_time.csv",FileType="spreadsheet")';
tumor = readmatrix("tumor.csv",FileType="spreadsheet")';
e25t1 = readmatrix("e25t1.csv",FileType="spreadsheet")';
e6t1 = readmatrix("e6t1.csv",FileType="spreadsheet")';
e1t1 = readmatrix("e1t1.csv",FileType="spreadsheet")';
patients = readcell("patCSV.csv"); patients = cell2table(patients);
cytolysis = readmatrix("cytolysis.csv",FileType="spreadsheet")';
e25 = readmatrix("e25.csv",FileType="spreadsheet")';
e6 = readmatrix("e6.csv",FileType="spreadsheet")';
e1 = readmatrix("e1.csv",FileType="spreadsheet")';
Ein = readmatrix("Ein.csv",FileType="spreadsheet")';
e25t1_sim = readmatrix("e25t1_sim_CI.csv",FileType="spreadsheet")';
e6t1_sim = readmatrix("e6t1_sim_CI.csv",FileType="spreadsheet")';
e1t1_sim = readmatrix("e1t1_sim_CI.csv",FileType="spreadsheet")';
tumor_sim = readmatrix("target_sim_CI.csv",FileType="spreadsheet")';

%% Append and filter datasets

s_var = cell2table(s_results(2:end,1:11),"VariableNames",s_results(1,:));
s_var = renamevars(s_var,"b","invb");
s_var.invb= 1./s_var.invb;
s_var.sl = s_var.s.^s_var.l;
p_var = cell2table(p_data(2:end,1:77),"VariableNames",p_data(1,1:77));
[var,il,ir] = innerjoin(s_var,p_var,"LeftKeys","Folder","RightKeys","Patient");
var = renamevars(var,"Folder","Patient");

disp("Ommited patients no patient data: " + newline); 
disp("p4" + newline); 
disp("p17" + newline);
disp("Ommited patients no xcelligence data: " + newline); 
disp("p2" + newline);

    % p2 (no xcelligence data), p4 (no patient data), p17 (no patient data)
    % are excluded

% exlude invalid fit datasets

SSR = var.loss;
high_SSR = find(SSR > 25);
high_SSR_p = var(high_SSR,"Patient");
disp("Ommited patients high SSR: " + newline); disp(high_SSR_p);
var(high_SSR,:) = [];

included_patients = cell2table(var.Patient);
[patients,ilp,irp] = innerjoin(patients,included_patients,"LeftKeys","patients","RightKeys","Var1");
time = time(:,ilp);
tumor = tumor(:,ilp);
e25t1 = e25t1(:,ilp);
e6t1 = e6t1(:,ilp);
e1t1 = e1t1(:,ilp);
cytolysis = cytolysis(:,ilp);
e25 = e25(:,ilp);
e6 = e6(:,ilp);
e1 = e1(:,ilp);
Ein = Ein(:,ilp);
e25t1_sim = e25t1_sim(:,ilp);
e6t1_sim = e6t1_sim(:,ilp);
e1t1_sim = e1t1_sim(:,ilp);
tumor_sim = tumor_sim(:,ilp);

% exclude invalid datasets

CI_difference = tumor(60,:) - tumor(12,:);
invalid_data = find(CI_difference < 0.0);
invalid_data_p = patients(invalid_data,"patients");
disp("Ommited patients invalid data: " + newline); disp(invalid_data_p);
var(invalid_data,:) = [];

included_patients = cell2table(var.Patient);
[patients,ilp,irp] = innerjoin(patients,included_patients,"LeftKeys","patients","RightKeys","Var1");
time = time(:,ilp);
tumor = tumor(:,ilp);
e25t1 = e25t1(:,ilp);
e6t1 = e6t1(:,ilp);
e1t1 = e1t1(:,ilp);
cytolysis = cytolysis(:,ilp);
e25 = e25(:,ilp);
e6 = e6(:,ilp);
e1 = e1(:,ilp);
Ein = Ein(:,ilp);
e25t1_sim = e25t1_sim(:,ilp);
e6t1_sim = e6t1_sim(:,ilp);
e1t1_sim = e1t1_sim(:,ilp);
tumor_sim = tumor_sim(:,ilp);

%exclude invalid datasets

CI_difference = tumor(12,:) - tumor(6,:);
invalid_data = find(CI_difference < 0.1);
invalid_data_p = patients(invalid_data,"patients");
disp("Ommited patients invalid data: " + newline); disp(invalid_data_p);
var(invalid_data,:) = [];

included_patients = cell2table(var.Patient);
[patients,ilp,irp] = innerjoin(patients,included_patients,"LeftKeys","patients","RightKeys","Var1");
time = time(:,ilp);
tumor = tumor(:,ilp);
e25t1 = e25t1(:,ilp);
e6t1 = e6t1(:,ilp);
e1t1 = e1t1(:,ilp);
cytolysis = cytolysis(:,ilp);
e25 = e25(:,ilp);
e6 = e6(:,ilp);
e1 = e1(:,ilp);
Ein = Ein(:,ilp);
e25t1_sim = e25t1_sim(:,ilp);
e6t1_sim = e6t1_sim(:,ilp);
e1t1_sim = e1t1_sim(:,ilp);
tumor_sim = tumor_sim(:,ilp);

%DLCBL

included_patients = var(strcmp(var.Disease, 'DLBCL'),'Patient');
[patients,ilp,irp] = innerjoin(patients,included_patients,"LeftKeys","patients","RightKeys","Patient");
time = time(:,ilp);
tumor = tumor(:,ilp);
e25t1 = e25t1(:,ilp);
e6t1 = e6t1(:,ilp);
e1t1 = e1t1(:,ilp);
cytolysis = cytolysis(:,ilp);
e25 = e25(:,ilp);
e6 = e6(:,ilp);
e1 = e1(:,ilp);
Ein = Ein(:,ilp);
e25t1_sim = e25t1_sim(:,ilp);
e6t1_sim = e6t1_sim(:,ilp);
e1t1_sim = e1t1_sim(:,ilp);
tumor_sim = tumor_sim(:,ilp);
var = var(strcmp(var.Disease, 'DLBCL'),:);

% means

[std_time,mean_time] = std(time,0,2,"omitnan");
[std_tumor,mean_tumor] = std(tumor,0,2,"omitnan");
[std_e25t1,mean_e25t1] = std(e25t1,0,2,"omitnan");
[std_e6t1,mean_e6t1] = std(e6t1,0,2,"omitnan");
[std_e1t1,mean_e1t1] = std(e1t1,0,2,"omitnan");
[std_cytolysis,mean_cytolysis] = std(cytolysis,0,2,"omitnan");
[std_e25,mean_e25] = std(e25,0,2,"omitnan");
[std_e6,mean_e6] = std(e6,0,2,"omitnan");
[std_e1,mean_e1] = std(e1,0,2,"omitnan");
[std_e25t1_sim,mean_e25t1_sim] = std(e25t1_sim,0,2,"omitnan");
[std_e6t1_sim,mean_e6t1_sim] = std(e6t1_sim,0,2,"omitnan");
[std_e1t1_sim,mean_e1t1_sim] = std(e1t1_sim,0,2,"omitnan");

% save datasets

save("mean_time.mat","mean_time");
save("mean_tumor.mat","mean_tumor");
save("mean_e25t1.mat","mean_e25t1");
save("mean_e6t1.mat","mean_e6t1");
save("mean_e1t1.mat","mean_e1t1");
save("mean_cytolysis.mat","mean_cytolysis");
save("mean_e25.mat","mean_e25");
save("mean_e6.mat","mean_e6");
save("mean_e1.mat","mean_e1");
save("std_tumor.mat","std_tumor");
save("std_e25t1.mat","std_e25t1");
save("std_e6t1.mat","std_e6t1");
save("std_e1t1.mat","std_e1t1");
save("std_cytolysis.mat","std_cytolysis");
save("std_e25.mat","std_e25");
save("std_e6.mat","std_e6");
save("std_e1.mat","std_e1");

save("time.mat","time");
save("tumor.mat","tumor");
save("e25t1.mat","e25t1");
save("e6t1.mat","e6t1");
save("e1t1.mat","e1t1");
save("cytolysis.mat","cytolysis");
save("e25.mat","e25");
save("e6.mat","e6");
save("e1.mat","e1");

save("mean_e25t1_sim.mat","mean_e25t1_sim");
save("mean_e6t1_sim.mat","mean_e6t1_sim");
save("mean_e1t1_sim.mat","mean_e1t1_sim");
save("e25t1_sim.mat","e25t1_sim");
save("e6t1_sim.mat","e6t1_sim");
save("e1t1_sim.mat","e1t1_sim");
save("tumor_sim.mat","tumor_sim");

% 28 day response category

CR_idx = find(strcmp(var.Day28_response,'CR') == 1);
PR_idx = find(strcmp(var.Day28_response,'PR') == 1);
PD_idx = find(strcmp(var.Day28_response,'PD') == 1);
Day28_response2(CR_idx) = "Response";
Day28_response2(PR_idx) = "Response";
Day28_response2(PD_idx) = "No response";
var = addvars(var,Day28_response2','NewVariableNames','Day28_response2');

% 90 day response category

CR_idx2 = find(strcmp(var.Day90_response,'CR') == 1);
PR_idx2 = find(strcmp(var.Day90_response,'PR') == 1);
PD_idx2 = find(strcmp(var.Day90_response,'PD') == 1);
ND_idx2 = find(strcmp(var.Day90_response,'ND') == 1);
Day90_response2(CR_idx2) = "Response";
Day90_response2(PR_idx2) = "Response";
Day90_response2(PD_idx2) = "No response";
Day90_response2(ND_idx2) = "No response";
var = addvars(var,Day90_response2','NewVariableNames','Day90_response2');

%% Figures

f = 0;

parameters = s_var.Properties.VariableNames(3:11);
parameters_order(1) = parameters(1);
parameters_order(2) = parameters(2);
parameters_order(3) = parameters(5);
parameters_order(4) = parameters(9);
parameters_order(5) = parameters(8);
parameters_order(6) = parameters(6);
parameters_order(7) = parameters(4);
parameters_order(8) = parameters(3);
parameters_order(9) = parameters(7);
parameters = parameters_order;

par_yscale(1) = "linear";
par_yscale(2) = "log";
par_yscale(3) = "linear";
par_yscale(4) = "linear";
par_yscale(5) = "linear";
par_yscale(6) = "linear";
par_yscale(7) = "log";
par_yscale(8) = "log";
par_yscale(9) = "linear";
par_fontsize = 24;
iqrmul = 100;
m = 1;
n = size(parameters,2);
par_plot_width = 18;
par_plot_height = 4.5;
ylim_a = [-0.05 0.37]; ylim_invb = [1e4 1e6]; ylim_kc = [-0.05 0.37];
ylim_s = [-0.90 4]; ylim_l = [-0.90 4]; ylim_kp2 = [-0.05 0.37];
ylim_k = [1e2 1e6]; ylim_kappa = [1e4 1e6]; ylim_kd = [-0.005 0.03]; ylim_sl = [0 10];
ylim_par = [ylim_a; ylim_invb; ylim_kc; ylim_s; ylim_l; ylim_kp2; ylim_k; ylim_kappa; ylim_kd];
y_exp_par = [0 0 0 0 0 0 0 0 0];
x_tick_angle_par = [45 0 0 0 0 0 0 0 0 ];
par_names = {'k_{p1}','C_T','k_c','K_{mr}','n_{ }','k_{p2}','K_{mp}','C_E','k_{d}'};

% %% All data plots
% 
% % Create data plots
% 
% xlim_time = [0 72];
% ylim_CI = [0 4.1];
% xticks_time = [0 24 48 72];
% yticks_CI = [0 1 2 3 4];
% plot_fontsize = 16;
% plot_fontsize2 = 22;
% legend_fontsize = 20;
% 
% % % % all datasets
% % % 
% % % f = f + 1;
% % % figure(f)
% % % set(gcf,'Units','inches','Position',[0.5,0.5,16,4])
% % % startColor = [0.5, 0.5, 0.5];
% % % 
% % % subplot(1,4,1)
% % % plot(mean_time,mean_tumor,'--k','Linewidth',2.5); hold on;
% % % patch([mean_time' flip(mean_time')],...
% % %       [(mean_tumor + std_tumor)' flip((mean_tumor - std_tumor)')],...
% % %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % % for ii = 1:size(time,2)
% % %     plot(time(:,ii),tumor(:,ii),'Color', startColor + ii/100);
% % % end
% % % xlim(xlim_time);
% % % ylim(ylim_CI);
% % % xticks(xticks_time);
% % % xlabel('Time (hours)'); ylabel('CI');
% % % set(gca,'FontSize',20)
% % % title("tumor")
% % % 
% % % subplot(1,4,4)
% % % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % % patch([mean_time' flip(mean_time')],...
% % %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% % %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % % for ii = 1:size(time,2)
% % %     plot(time(:,ii),e25t1(:,ii),'Color', startColor + ii/100);
% % % end
% % % xlim(xlim_time);
% % % ylim(ylim_CI);
% % % xticks(xticks_time);
% % % xlabel('Time (hours)'); ylabel('CI');
% % % set(gca,'FontSize',20)
% % % title("E:T 25:1")
% % % 
% % % subplot(1,4,3)
% % % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % % patch([mean_time' flip(mean_time')],...
% % %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% % %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % % for ii = 1:size(time,2)
% % %     plot(time(:,ii),e6t1(:,ii),'Color', startColor + ii/100);
% % % end
% % % xlim(xlim_time);
% % % ylim(ylim_CI);
% % % xticks(xticks_time);
% % % xlabel('Time (hours)'); ylabel('CI');
% % % set(gca,'FontSize',20)
% % % title("E:T 6.25:1")
% % % 
% % % subplot(1,4,2)
% % % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % % patch([mean_time' flip(mean_time')],...
% % %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% % %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % % for ii = 1:size(time,2)
% % %     plot(time(:,ii),e1t1(:,ii),'Color', startColor + ii/100);
% % % end
% % % xlim(xlim_time);
% % % ylim(ylim_CI);
% % % xticks(xticks_time);
% % % xlabel('Time (hours)'); ylabel('CI');
% % % set(gca,'FontSize',20)
% % % title("E:T 1:1")
% % % 
% % % exportgraphics(gcf,'all_data.tiff','Resolution',1200)
% % 
% % 28 day response
% 
% % data CR
% patients_CR = var(strcmp(var.Day28_response2,'Response'),"Patient");
% [patients_CR,ilpcr,irpcr] = innerjoin(patients,patients_CR,"LeftKeys","patients","RightKeys","Patient");
% time_CR = time(:,ilpcr);
% tumor_CR = tumor(:,ilpcr);
% e25t1_CR = e25t1(:,ilpcr);
% e6t1_CR = e6t1(:,ilpcr);
% e1t1_CR = e1t1(:,ilpcr);
% 
% % means CR
% [std_time_CR,mean_time_CR] = std(time_CR,0,2,"omitnan");
% [std_tumor_CR,mean_tumor_CR] = std(tumor_CR,0,2,"omitnan");
% [std_e25t1_CR,mean_e25t1_CR] = std(e25t1_CR,0,2,"omitnan");
% [std_e6t1_CR,mean_e6t1_CR] = std(e6t1_CR,0,2,"omitnan");
% [std_e1t1_CR,mean_e1t1_CR] = std(e1t1_CR,0,2,"omitnan");
% 
% % data NCR
% patients_NCR = var(strcmp(var.Day28_response2,'No response'),"Patient");
% [patients_NCR,ilpncr,irpncr] = innerjoin(patients,patients_NCR,"LeftKeys","patients","RightKeys","Patient");
% time_NCR = time(:,ilpncr);
% tumor_NCR = tumor(:,ilpncr);
% e25t1_NCR = e25t1(:,ilpncr);
% e6t1_NCR = e6t1(:,ilpncr);
% e1t1_NCR = e1t1(:,ilpncr);
% 
% % means NCR
% [std_time_NCR,mean_time_NCR] = std(time_NCR,0,2,"omitnan");
% [std_tumor_NCR,mean_tumor_NCR] = std(tumor_NCR,0,2,"omitnan");
% [std_e25t1_NCR,mean_e25t1_NCR] = std(e25t1_NCR,0,2,"omitnan");
% [std_e6t1_NCR,mean_e6t1_NCR] = std(e6t1_NCR,0,2,"omitnan");
% [std_e1t1_NCR,mean_e1t1_NCR] = std(e1t1_NCR,0,2,"omitnan");
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,12,12])
% startColor = [0.5, 0.5, 0.5];
% %sgtitle("cytotoxicity by day 28 response (blue - CR, red - NCR)")
% m = 3;
% n = 3;
% 
% subplot(m,n,1)
% % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e1t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e1t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 1:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(m,n,2)
% % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e6t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e6t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 6.25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(m,n,3)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     x1 = (categorical_data == 'Response'); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); par_n_x2 = sum(x2);
% 
% % subplot(m,n,4)
% % % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % % patch([mean_time' flip(mean_time')],...
% % %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% % %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % xlim(xlim_time);
% % ylim(ylim_CI);
% % xticks(xticks_time);
% % yticks(yticks_CI);
% % axis([10,11,10,11])
% % H = legend( ['Response (R) (n=',num2str(par_n_x1),')'],...
% %             'SD R',...
% %             ['No response (NR) (n=',num2str(par_n_x2),')'],...
% %             'SD NR',...
% %             'Location','west');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% % axis off
% 
% %each time course
% 
% startColor = [0.0, 0.0, 0.5];
% 
% subplot(m,n,6)
% 
% plot(mean_time_CR,mean_e25t1_CR,'--b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_CR,2)
%     plot(time_CR(:,ii),e25t1_CR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,5)
% 
% plot(mean_time_CR,mean_e6t1_CR,'--b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_CR,2)
%     plot(time_CR(:,ii),e6t1_CR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,4)
% 
% plot(mean_time_CR,mean_e1t1_CR,'--b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_CR,2)
%     plot(time_CR(:,ii),e1t1_CR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% startColor = [0.5, 0.0, 0.0];
% 
% subplot(m,n,9)
% 
% plot(mean_time_NCR,mean_e25t1_NCR,'--r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_NCR,2)
%     plot(time_NCR(:,ii),e25t1_NCR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,8)
% 
% plot(mean_time_NCR,mean_e6t1_NCR,'--r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_NCR,2)
%     plot(time_NCR(:,ii),e6t1_NCR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,7)
% 
% plot(mean_time_NCR,mean_e1t1_NCR,'--r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_NCR,2)
%     plot(time_NCR(:,ii),e1t1_NCR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% exportgraphics(gcf,'all_data_day28response_tc.tiff','Resolution',1200)

% % 90 day response
% 
% % data CR
% patients_CR = var(strcmp(var.Day90_response2,'Response'),"Patient");
% [patients_CR,ilpcr,irpcr] = innerjoin(patients,patients_CR,"LeftKeys","patients","RightKeys","Patient");
% time_CR = time(:,ilpcr);
% tumor_CR = tumor(:,ilpcr);
% e25t1_CR = e25t1(:,ilpcr);
% e6t1_CR = e6t1(:,ilpcr);
% e1t1_CR = e1t1(:,ilpcr);
% 
% % means CR
% [std_time_CR,mean_time_CR] = std(time_CR,0,2,"omitnan");
% [std_tumor_CR,mean_tumor_CR] = std(tumor_CR,0,2,"omitnan");
% [std_e25t1_CR,mean_e25t1_CR] = std(e25t1_CR,0,2,"omitnan");
% [std_e6t1_CR,mean_e6t1_CR] = std(e6t1_CR,0,2,"omitnan");
% [std_e1t1_CR,mean_e1t1_CR] = std(e1t1_CR,0,2,"omitnan");
% 
% % data NCR
% patients_NCR = var(strcmp(var.Day90_response2,'No response'),"Patient");
% [patients_NCR,ilpncr,irpncr] = innerjoin(patients,patients_NCR,"LeftKeys","patients","RightKeys","Patient");
% time_NCR = time(:,ilpncr);
% tumor_NCR = tumor(:,ilpncr);
% e25t1_NCR = e25t1(:,ilpncr);
% e6t1_NCR = e6t1(:,ilpncr);
% e1t1_NCR = e1t1(:,ilpncr);
% 
% % means NCR
% [std_time_NCR,mean_time_NCR] = std(time_NCR,0,2,"omitnan");
% [std_tumor_NCR,mean_tumor_NCR] = std(tumor_NCR,0,2,"omitnan");
% [std_e25t1_NCR,mean_e25t1_NCR] = std(e25t1_NCR,0,2,"omitnan");
% [std_e6t1_NCR,mean_e6t1_NCR] = std(e6t1_NCR,0,2,"omitnan");
% [std_e1t1_NCR,mean_e1t1_NCR] = std(e1t1_NCR,0,2,"omitnan");
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,16,4])
% startColor = [0.5, 0.5, 0.5];
% %sgtitle("cytotoxicity by day 28 response (blue - CR, red - NCR)")
% 
% subplot(1,4,1)
% % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e1t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e1t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 1:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,2)
% % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e6t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e6t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 6.25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,3)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
%     categorical_data = categorical(table2array(var(:,'Day90_response2')));
%     x1 = (categorical_data == 'Response'); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); par_n_x2 = sum(x2);
% 
% subplot(1,4,4)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% axis([10,11,10,11])
% H = legend( ['Response (R) (n=',num2str(par_n_x1),')'],...
%             'SD R',...
%             ['No response (NR) (n=',num2str(par_n_x2),')'],...
%             'SD NR',...
%             'Location','west');
% set(H,'FontSize',legend_fontsize)
% legend boxoff
% axis off
% 
% % each time course
% 
% % subplot(3,3,4)
% % 
% % plot(mean_time_CR,mean_e25t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e25t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,7)
% % 
% % plot(mean_time_NCR,mean_e25t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e25t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,5)
% % 
% % plot(mean_time_CR,mean_e6t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e6t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,8)
% % 
% % plot(mean_time_NCR,mean_e6t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e6t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,6)
% % 
% % plot(mean_time_CR,mean_e1t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e1t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,9)
% % 
% % plot(mean_time_NCR,mean_e1t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e1t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% 
% exportgraphics(gcf,'all_data_day90response.tiff','Resolution',1200)
% 
% % Relapse
% 
% % data CR
% patients_CR = var(strcmp(var.Relapse,'No'),"Patient");
% [patients_CR,ilpcr,irpcr] = innerjoin(patients,patients_CR,"LeftKeys","patients","RightKeys","Patient");
% time_CR = time(:,ilpcr);
% tumor_CR = tumor(:,ilpcr);
% e25t1_CR = e25t1(:,ilpcr);
% e6t1_CR = e6t1(:,ilpcr);
% e1t1_CR = e1t1(:,ilpcr);
% 
% % means CR
% [std_time_CR,mean_time_CR] = std(time_CR,0,2,"omitnan");
% [std_tumor_CR,mean_tumor_CR] = std(tumor_CR,0,2,"omitnan");
% [std_e25t1_CR,mean_e25t1_CR] = std(e25t1_CR,0,2,"omitnan");
% [std_e6t1_CR,mean_e6t1_CR] = std(e6t1_CR,0,2,"omitnan");
% [std_e1t1_CR,mean_e1t1_CR] = std(e1t1_CR,0,2,"omitnan");
% 
% % data NCR
% patients_NCR = var(strcmp(var.Relapse,'Yes'),"Patient");
% [patients_NCR,ilpncr,irpncr] = innerjoin(patients,patients_NCR,"LeftKeys","patients","RightKeys","Patient");
% time_NCR = time(:,ilpncr);
% tumor_NCR = tumor(:,ilpncr);
% e25t1_NCR = e25t1(:,ilpncr);
% e6t1_NCR = e6t1(:,ilpncr);
% e1t1_NCR = e1t1(:,ilpncr);
% 
% % means NCR
% [std_time_NCR,mean_time_NCR] = std(time_NCR,0,2,"omitnan");
% [std_tumor_NCR,mean_tumor_NCR] = std(tumor_NCR,0,2,"omitnan");
% [std_e25t1_NCR,mean_e25t1_NCR] = std(e25t1_NCR,0,2,"omitnan");
% [std_e6t1_NCR,mean_e6t1_NCR] = std(e6t1_NCR,0,2,"omitnan");
% [std_e1t1_NCR,mean_e1t1_NCR] = std(e1t1_NCR,0,2,"omitnan");
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,16,4])
% startColor = [0.5, 0.5, 0.5];
% %sgtitle("cytotoxicity by day 28 response (blue - CR, red - NCR)")
% 
% subplot(1,4,1)
% % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e1t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e1t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 1:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,2)
% % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e6t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e6t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 6.25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,3)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     x1 = (categorical_data == 'No'); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); par_n_x2 = sum(x2);
% 
% subplot(1,4,4)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% axis([10,11,10,11])
% H = legend( ['No Relapse (N) (n=',num2str(par_n_x1),')'],...
%             'SD N',...
%             ['Relapse (Y) (n=',num2str(par_n_x2),')'],...
%             'SD Y',...
%             'Location','west');
% set(H,'FontSize',legend_fontsize)
% legend boxoff
% axis off
% 
% % each time course
% 
% % subplot(3,3,4)
% % 
% % plot(mean_time_CR,mean_e25t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e25t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,7)
% % 
% % plot(mean_time_NCR,mean_e25t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e25t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,5)
% % 
% % plot(mean_time_CR,mean_e6t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e6t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,8)
% % 
% % plot(mean_time_NCR,mean_e6t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e6t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,6)
% % 
% % plot(mean_time_CR,mean_e1t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e1t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,9)
% % 
% % plot(mean_time_NCR,mean_e1t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e1t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% 
% exportgraphics(gcf,'all_data_relapse.tiff','Resolution',1200)
% 
%
% %disease
% 
% % data DLBCL
% patients_DLBCL = var(strcmp(var.Disease,'DLBCL'),"Patient");
% [patients_DLBCL,ilpdlbcl,irpdlbcl] = innerjoin(patients,patients_DLBCL,"LeftKeys","patients","RightKeys","Patient");
% time_DLBCL = time(:,ilpdlbcl);
% tumor_DLBCL = tumor(:,ilpdlbcl);
% e25t1_DLBCL = e25t1(:,ilpdlbcl);
% e6t1_DLBCL = e6t1(:,ilpdlbcl);
% e1t1_DLBCL = e1t1(:,ilpdlbcl);
% 
% % means DLBCL
% [std_time_DLBCL,mean_time_DLBCL] = std(time_DLBCL,0,2,"omitnan");
% [std_tumor_DLBCL,mean_tumor_DLBCL] = std(tumor_DLBCL,0,2,"omitnan");
% [std_e25t1_DLBCL,mean_e25t1_DLBCL] = std(e25t1_DLBCL,0,2,"omitnan");
% [std_e6t1_DLBCL,mean_e6t1_DLBCL] = std(e6t1_DLBCL,0,2,"omitnan");
% [std_e1t1_DLBCL,mean_e1t1_DLBCL] = std(e1t1_DLBCL,0,2,"omitnan");
% 
% % data MCL
% patients_MCL = var(strcmp(var.Disease,'MCL'),"Patient");
% [patients_MCL,ilpMCL,irpMCL] = innerjoin(patients,patients_MCL,"LeftKeys","patients","RightKeys","Patient");
% time_MCL = time(:,ilpMCL);
% tumor_MCL = tumor(:,ilpMCL);
% e25t1_MCL = e25t1(:,ilpMCL);
% e6t1_MCL = e6t1(:,ilpMCL);
% e1t1_MCL = e1t1(:,ilpMCL);
% 
% % means MCL
% [std_time_MCL,mean_time_MCL] = std(time_MCL,0,2,"omitnan");
% [std_tumor_MCL,mean_tumor_MCL] = std(tumor_MCL,0,2,"omitnan");
% [std_e25t1_MCL,mean_e25t1_MCL] = std(e25t1_MCL,0,2,"omitnan");
% [std_e6t1_MCL,mean_e6t1_MCL] = std(e6t1_MCL,0,2,"omitnan");
% [std_e1t1_MCL,mean_e1t1_MCL] = std(e1t1_MCL,0,2,"omitnan");
% 
% % data FL
% patients_FL = var(strcmp(var.Disease,'FL'),"Patient");
% [patients_FL,ilpFL,irpFL] = innerjoin(patients,patients_FL,"LeftKeys","patients","RightKeys","Patient");
% time_FL = time(:,ilpFL);
% tumor_FL = tumor(:,ilpFL);
% e25t1_FL = e25t1(:,ilpFL);
% e6t1_FL = e6t1(:,ilpFL);
% e1t1_FL = e1t1(:,ilpFL);
% 
% % means FL
% [std_time_FL,mean_time_FL] = std(time_FL,0,2,"omitnan");
% [std_tumor_FL,mean_tumor_FL] = std(tumor_FL,0,2,"omitnan");
% [std_e25t1_FL,mean_e25t1_FL] = std(e25t1_FL,0,2,"omitnan");
% [std_e6t1_FL,mean_e6t1_FL] = std(e6t1_FL,0,2,"omitnan");
% [std_e1t1_FL,mean_e1t1_FL] = std(e1t1_FL,0,2,"omitnan");
% 
% % data CLL
% patients_CLL = var(strcmp(var.Disease,'CLL'),"Patient");
% [patients_CLL,ilpCLL,irpCLL] = innerjoin(patients,patients_CLL,"LeftKeys","patients","RightKeys","Patient");
% time_CLL = time(:,ilpCLL);
% tumor_CLL = tumor(:,ilpCLL);
% e25t1_CLL = e25t1(:,ilpCLL);
% e6t1_CLL = e6t1(:,ilpCLL);
% e1t1_CLL = e1t1(:,ilpCLL);
% 
% % means CLL
% [std_time_CLL,mean_time_CLL] = std(time_CLL,0,2,"omitnan");
% [std_tumor_CLL,mean_tumor_CLL] = std(tumor_CLL,0,2,"omitnan");
% [std_e25t1_CLL,mean_e25t1_CLL] = std(e25t1_CLL,0,2,"omitnan");
% [std_e6t1_CLL,mean_e6t1_CLL] = std(e6t1_CLL,0,2,"omitnan");
% [std_e1t1_CLL,mean_e1t1_CLL] = std(e1t1_CLL,0,2,"omitnan");
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,16,4])
% startColor = [0.5, 0.5, 0.5];
% %sgtitle("cytotoxicity by disease")
% 
% subplot(1,4,1)
% % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_DLBCL,mean_e1t1_DLBCL,'-g','Linewidth',2.5); hold on;
% patch([mean_time_DLBCL' flip(mean_time_DLBCL')],...
%       [(mean_e1t1_DLBCL + std_e1t1_DLBCL)' flip((mean_e1t1_DLBCL - std_e1t1_DLBCL)')],...
%       [0 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_MCL,mean_e1t1_MCL,'-y','Linewidth',2.5,'Color',[0.9290 0.6940 0.1250]); hold on;
% patch([mean_time_MCL' flip(mean_time_MCL')],...
%       [(mean_e1t1_MCL + std_e1t1_MCL)' flip((mean_e1t1_MCL - std_e1t1_MCL)')],...
%       [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_FL,mean_e1t1_FL,'-c','Linewidth',2.5); hold on;
% patch([mean_time_FL' flip(mean_time_FL')],...
%       [(mean_e1t1_FL + std_e1t1_FL)' flip((mean_e1t1_FL - std_e1t1_FL)')],...
%       [0 1 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CLL,mean_e1t1_CLL,'-m','Linewidth',2.5); hold on;
% patch([mean_time_CLL' flip(mean_time_CLL')],...
%       [(mean_e1t1_CLL + std_e1t1_CLL)' flip((mean_e1t1_CLL - std_e1t1_CLL)')],...
%       [1 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 1:1")
% % H = legend('all','SD all', 'DLBCL',...
% %             'SD DLBCL','MCL', 'SD MCL',...
% %             'FL', 'SD FL', 'CLL', 'SD CLL',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,2)
% % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_DLBCL,mean_e6t1_DLBCL,'-g','Linewidth',2.5); hold on;
% patch([mean_time_DLBCL' flip(mean_time_DLBCL')],...
%       [(mean_e6t1_DLBCL + std_e6t1_DLBCL)' flip((mean_e6t1_DLBCL - std_e6t1_DLBCL)')],...
%       [0 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_MCL,mean_e6t1_MCL,'-y','Linewidth',2.5,'Color',[0.9290 0.6940 0.1250]); hold on;
% patch([mean_time_MCL' flip(mean_time_MCL')],...
%       [(mean_e6t1_MCL + std_e6t1_MCL)' flip((mean_e6t1_MCL - std_e6t1_MCL)')],...
%       [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_FL,mean_e6t1_FL,'-c','Linewidth',2.5); hold on;
% patch([mean_time_FL' flip(mean_time_FL')],...
%       [(mean_e6t1_FL + std_e6t1_FL)' flip((mean_e6t1_FL - std_e6t1_FL)')],...
%       [0 1 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CLL,mean_e6t1_CLL,'-m','Linewidth',2.5); hold on;
% patch([mean_time_CLL' flip(mean_time_CLL')],...
%       [(mean_e6t1_CLL + std_e6t1_CLL)' flip((mean_e6t1_CLL - std_e6t1_CLL)')],...
%       [1 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 6.25:1")
% % H = legend('all','SD all', 'DLBCL',...
% %             'SD DLBCL','MCL', 'SD MCL',...
% %             'FL', 'SD FL', 'CLL', 'SD CLL',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,3)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_DLBCL,mean_e25t1_DLBCL,'-g','Linewidth',2.5); hold on;
% patch([mean_time_DLBCL' flip(mean_time_DLBCL')],...
%       [(mean_e25t1_DLBCL + std_e25t1_DLBCL)' flip((mean_e25t1_DLBCL - std_e25t1_DLBCL)')],...
%       [0 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_MCL,mean_e25t1_MCL,'-y','Linewidth',2.5,'Color',[0.9290 0.6940 0.1250]); hold on;
% patch([mean_time_MCL' flip(mean_time_MCL')],...
%       [(mean_e25t1_MCL + std_e25t1_MCL)' flip((mean_e25t1_MCL - std_e25t1_MCL)')],...
%       [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_FL,mean_e25t1_FL,'-c','Linewidth',2.5); hold on;
% patch([mean_time_FL' flip(mean_time_FL')],...
%       [(mean_e25t1_FL + std_e25t1_FL)' flip((mean_e25t1_FL - std_e25t1_FL)')],...
%       [0 1 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CLL,mean_e25t1_CLL,'-m','Linewidth',2.5); hold on;
% patch([mean_time_CLL' flip(mean_time_CLL')],...
%       [(mean_e25t1_CLL + std_e25t1_CLL)' flip((mean_e25t1_CLL - std_e25t1_CLL)')],...
%       [1 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 25:1")
% % H = legend('all','SD all', 'DLBCL',...
% %             'SD DLBCL','MCL', 'SD MCL',...
% %             'FL', 'SD FL', 'CLL', 'SD CLL',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
%     categorical_data = categorical(table2array(var(:,'Disease')));
%     x1 = (categorical_data == 'DLBCL'); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'MCL'); par_n_x2 = sum(x2);
%     x3 = (categorical_data == 'FL'); par_n_x3 = sum(x3);
%     x4 = (categorical_data == 'CLL'); par_n_x4 = sum(x4);
% 
% subplot(1,4,4)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_DLBCL,mean_e25t1_DLBCL,'-g','Linewidth',2.5); hold on;
% patch([mean_time_DLBCL' flip(mean_time_DLBCL')],...
%       [(mean_e25t1_DLBCL + std_e25t1_DLBCL)' flip((mean_e25t1_DLBCL - std_e25t1_DLBCL)')],...
%       [0 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_MCL,mean_e25t1_MCL,'-y','Linewidth',2.5,'Color',[0.9290 0.6940 0.1250]); hold on;
% patch([mean_time_MCL' flip(mean_time_MCL')],...
%       [(mean_e25t1_MCL + std_e25t1_MCL)' flip((mean_e25t1_MCL - std_e25t1_MCL)')],...
%       [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_FL,mean_e25t1_FL,'-c','Linewidth',2.5); hold on;
% patch([mean_time_FL' flip(mean_time_FL')],...
%       [(mean_e25t1_FL + std_e25t1_FL)' flip((mean_e25t1_FL - std_e25t1_FL)')],...
%       [0 1 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CLL,mean_e25t1_CLL,'-m','Linewidth',2.5); hold on;
% patch([mean_time_CLL' flip(mean_time_CLL')],...
%       [(mean_e25t1_CLL + std_e25t1_CLL)' flip((mean_e25t1_CLL - std_e25t1_CLL)')],...
%       [1 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% axis([10,11,10,11])
% H = legend( ['DLBCL (n=',num2str(par_n_x1),')'],...
%             'SD DLBCL',...
%             ['MCL (n=',num2str(par_n_x2),')'],...
%             'SD MCL',...
%              ['FL (n=',num2str(par_n_x3),')'],...
%             'SD FL',...
%              ['CLL (n=',num2str(par_n_x4),')'],...
%             'SD CLL',...
%             'Location','west');
% set(H,'FontSize',legend_fontsize)
% legend boxoff
% axis off
% 
% exportgraphics(gcf,'all_data_disease.tiff','Resolution',1200)

%% Create parameter plots

% SSR bar
f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,4,4])

histogram(var.loss,10); hold on;
SSR_mean = mean(var.loss);
xline(SSR_mean,'--r','LineWidth',2);
ylabel("# datasets");
xlabel("SSR");
set(gca,'FontSize',par_fontsize)

% Parameter scatter
f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width,par_plot_height])
%sgtitle("Parameter Distribution")
t = tiledlayout(1,9);

for i = 1:size(parameters,2)

    nexttile
    parameter_data = table2array(var(:,parameters(i)));
    param = categorical(cellstr(par_names(i)));
    x = ones([size(parameter_data,1),1]);
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];

    boxplot(parameter_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=2, ...
            Whisker=3); 

    
    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));

    hold on;
    [par_std(i), par_mean(i)] = std(parameter_data);
    par_n = sum(x(:,1));
    plot(par_mean(i),'ok','MarkerSize',7,...
                      'MarkerEdgeColor','k',...
                      'MarkerFaceColor','k');
    swarmchart(x,parameter_data,[],[0.5 0.5 0.5]);
    
    if i == 1
        txt = strcat(string(param),' (n=', num2str(par_n),')');
        xlabel(txt)%,'Interpreter','none')
    else
        xlabel(param)%,'Interpreter','none')
    end
    
    hold off;

end

exportgraphics(gcf,'parameters_dist.tiff','Resolution',1200)

% Parameter scatter day 28 response
parameters2 = s_var.Properties.VariableNames(5:11);
parameters_order2(1) = parameters2(3);
parameters_order2(2) = parameters2(7);
parameters_order2(3) = parameters2(6);
parameters_order2(4) = parameters2(4);
parameters_order2(5) = parameters2(2);
parameters_order2(6) = parameters2(1);
parameters2 = parameters_order2;

par_yscale(1) = "linear";
par_yscale(2) = "linear";
par_yscale(3) = "linear";
par_yscale(4) = "linear";
par_yscale(5) = "log";
par_yscale(6) = "log";
ylim_par = [ylim_kc; ylim_s; ylim_l; ylim_kp2; ylim_k; ylim_kappa];
y_exp_par = [0 0 0 0 0 0];
x_tick_angle_par = [0 0 0 0 0 0];
par_names = {'k_c','K_{mr}','n_{ }','k_{p2}','K_{mp}','C_E'};
par_fontsize = 18;

f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/1.5,par_plot_height])
%sgtitle("Day 28 Response")
t = tiledlayout(1,6);

for i = 1:size(parameters2,2)

    %subplot(m,n,i)
    nexttile
    parameter_data = table2array(var(:,parameters2(i)));
    categorical_data = categorical(table2array(var(:,'Day28_response2')));
    param = categorical(cellstr(par_names(i)));
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];

    boxplot(parameter_data, categorical_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=0.5,...
            GroupOrder={'Response','No response'}); 

    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));

    
    hold on;
    x = zeros([size(categorical_data,1),1]);
    x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
    x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
    
    swarmchart(x(x1),parameter_data(x1),[],[0 0 1]);
    swarmchart(x(x2),parameter_data(x2),[],[1 0 0]);
    plot(1,par_mean_x1,'ob','MarkerSize',7,'MarkerFaceColor','b');
    plot(2,par_mean_x2,'or','MarkerSize',7,'MarkerFaceColor','r');

%     if i == 8
%         xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     else
%         xticklabels({'R','NR'})
%     end

    xlabel(param)%,'Interpreter','none')
    xticklabels({'R','NR'})
    xtickangle(0)

    % statistics for each group
    [chu(i,1),cpu(i,1),stats] = ranksum(parameter_data(x1),parameter_data(x2));
    stats_28d_CI(i,:) = meanEffectSize(parameter_data(x1),parameter_data(x2),...
                        "VarianceType","unequal",ConfidenceIntervalType="bootstrap");

    hold off;

end

exportgraphics(gcf,'parameters_dist_28dayresponse.tiff','Resolution',1200)

% Parameter scatter day 90 response

f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/1.5,par_plot_height])
%sgtitle("Day 90 Response")
t = tiledlayout(1,6);

for i = 1:size(parameters2,2)

    %subplot(m,n,i)
    nexttile
    parameter_data = table2array(var(:,parameters2(i)));
    categorical_data = categorical(table2array(var(:,'Day90_response2')));
    missing = find(categorical_data == 'missing');
    parameter_data(missing) = [];
    categorical_data(missing) = [];
    param = categorical(cellstr(par_names(i)));
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];


    boxplot(parameter_data, categorical_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=0.5,...
            GroupOrder={'Response','No response'}); 

    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));
    xlabel(param)%,'Interpreter','none')

    hold on;

    x = zeros([size(categorical_data,1),1]);
    x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
    x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
    
    swarmchart(x(x1),parameter_data(x1),[],[0 0 1]);
    swarmchart(x(x2),parameter_data(x2),[],[1 0 0]);
    plot(1,par_mean_x1,'ob','MarkerSize',7,'MarkerFaceColor','b');
    plot(2,par_mean_x2,'or','MarkerSize',7,'MarkerFaceColor','r');

    if i == 8
        xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
    else
        xticklabels({'R','NR'})
    end

    xtickangle(0)

    % statistics for each group
    [chu(i,2),cpu(i,2)] = ranksum(parameter_data(x1),parameter_data(x2));
    stats_90d_CI(i,:) = meanEffectSize(parameter_data(x1),parameter_data(x2),...
                        "VarianceType","unequal",ConfidenceIntervalType="bootstrap");


    hold off;

end

exportgraphics(gcf,'parameters_dist_90dayresponse.tiff','Resolution',1200)

% Parameter scatter CRS

f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/1.5,par_plot_height])
%sgtitle("CRS")
t = tiledlayout(1,6);

for i = 1:size(parameters2,2)

    %subplot(m,n,i)
    nexttile
    parameter_data = table2array(var(:,parameters2(i)));
    categorical_data = categorical(table2array(var(:,'CRS')));
    param = categorical(cellstr(par_names(i)));
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];

    boxplot(parameter_data, categorical_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=0.5,...
            GroupOrder={'Yes','No'}); 

    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));
    xlabel(param)%,'Interpreter','none')
    
    hold on;
    x = zeros([size(categorical_data,1),1]);
    x1 = (categorical_data == 'Yes'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
    x2 = (categorical_data == 'No'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
    
    swarmchart(x,parameter_data,[],[0.5 0.5 0.5]);
    plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
    plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');

    if i == 1
        xticklabels({['Y (n=',num2str(par_n_x1),')'],['N (n=',num2str(par_n_x2),')']})
    else
        xticklabels({'Y','N'})
    end

    % statistics for each group
    [chu(i,4),cpu(i,4)] = ranksum(parameter_data(x1),parameter_data(x2));

    hold off;

end

exportgraphics(gcf,'parameters_dist_CRS.tiff','Resolution',1200)



% Parameter scatter Neurotoxicity

f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/1.5,par_plot_height])
t = tiledlayout(1,6);

for i = 1:size(parameters2,2)

    %subplot(m,n,i)
    nexttile
    parameter_data = table2array(var(:,parameters2(i)));
    categorical_data = categorical(table2array(var(:,'Neurotoxicity')));
    param = categorical(cellstr(par_names(i)));
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];

    boxplot(parameter_data, categorical_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=0.5,...
            GroupOrder={'Yes','No'}); 

    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));
    xlabel(param)%,'Interpreter','none')
    
    hold on;
    x = zeros([size(categorical_data,1),1]);
    x1 = (categorical_data == 'Yes'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
    x2 = (categorical_data == 'No'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
    
    swarmchart(x,parameter_data,[],[0.5 0.5 0.5]);
    plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
    plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');

    if i == 1
        xticklabels({['Y (n=',num2str(par_n_x1),')'],['N (n=',num2str(par_n_x2),')']})
    else
        xticklabels({'Y','N'})
    end

    % statistics for each group
    [chu(i,5),cpu(i,5)] = ranksum(parameter_data(x1),parameter_data(x2));

    hold off;

end

exportgraphics(gcf,'parameters_dist_neurotoxicity.tiff','Resolution',1200)


% Parameter scatter Manufacturing Days

f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/1.5,par_plot_height])
%sgtitle("Manufacturing Days")
t = tiledlayout(1,6);

for i = 1:size(parameters2,2)

    %subplot(m,n,i)
    nexttile
    parameter_data = table2array(var(:,parameters2(i)));
    categorical_data = categorical(table2array(var(:,'Manufacturing protocol - Harvest Day')));
    param = categorical(cellstr(par_names(i)));
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];

    boxplot(parameter_data, categorical_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=0.5,...
            GroupOrder={'8','12'}); 

    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));
    xlabel(param)%,'Interpreter','none')

    hold on;
    x = zeros([size(categorical_data,1),1]);
    x1 = (categorical_data == '8'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
    x2 = (categorical_data == '12'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
    
    swarmchart(x,parameter_data,[],[0.5 0.5 0.5]);
    plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
    plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');

    if i == 1
        xticklabels({['8 (n=',num2str(par_n_x1),')'],['12 (n=',num2str(par_n_x2),')']})
    else
        xticklabels({'8','12'})
    end

    % statistics for each group
    [chu(i,6),cpu(i,6)] = ranksum(parameter_data(x1),parameter_data(x2));

    hold off;

end

exportgraphics(gcf,'parameters_dist_manufacturingdays.tiff','Resolution',1200)


% % Responders / non-responders CD4:CD8 ratio
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cd4cd8_data = table2array(var(:,'CD4_CD8_Ratio_CART'));
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd4cd8_data)*iqrmul,iqr(cd4cd8_data)*iqrmul];
% 
%     boxplot(cd4cd8_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('CD4:CD8')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd4cd8_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd4cd8_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4cd8_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xtickangle(0)
% 
%     % statistics for each group
%     [cd4cd8ht(1),cd4cd8pt(1)] = ttest2(cd4cd8_data(x1),cd4cd8_data(x2),'Vartype','unequal');
%     [cd4cd8pu(1),cd4cd8hu(1)] = ranksum(cd4cd8_data(x1),cd4cd8_data(x2));
% 
% exportgraphics(gcf,'response_day28_CD4CD8.tiff','Resolution',1200)
% 
% 
% % Responders / non-responders CD4:CD8 ratio
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cd4cd8_data = table2array(var(:,'CD4_CD8_Ratio_CART'));
%     categorical_data = categorical(table2array(var(:,'Day90_response2')));
%     DL = [-iqr(cd4cd8_data)*iqrmul,iqr(cd4cd8_data)*iqrmul];
% 
%     boxplot(cd4cd8_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('CD4:CD8')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd4cd8_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd4cd8_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4cd8_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xtickangle(0)
% 
%     % statistics for each group
%     [cd4cd8ht(2),cd4cd8pt(2)] = ttest2(cd4cd8_data(x1),cd4cd8_data(x2),'Vartype','unequal');
%     [cd4cd8pu(2),cd4cd8hu(2)] = ranksum(cd4cd8_data(x1),cd4cd8_data(x2));
% 
% 
% exportgraphics(gcf,'response_day90_CD4CD8.tiff','Resolution',1200)
% 
% % Parameter scatter CD4:CD8 ratio
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width*0.6,par_plot_height*2])
% %sgtitle("CD4:CD8 ratio")
% m = 1;
% n = 7;
% t = tiledlayout(2,3);
% 
% 
% for i = 1:size(parameters2,2)
% 
%     %subplot(m,n,i)
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     cd4cd8_data = table2array(var(:,'CD4_CD8_Ratio_CART'));
%     param = categorical(cellstr(par_names(i)));
% 
%     scatter(cd4cd8_data, parameter_data,72,'kx','LineWidth',1.5);
% 
%     ylim(ylim_par(i,:));
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)  
%     xlim([0 22])
%     xlabel('CD4:CD8')
%     ylabel(param)%,'Interpreter','none')
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_CD4CD8.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4 Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd4tscm_data)*iqrmul,iqr(cd4tscm_data)*iqrmul];
% 
%     boxplot(cd4tscm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 Tscm')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd4tscm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd4tscm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4tscm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [tscmht(1),tscmpt(1)] = ttest2(cd4tscm_data(x1),cd4tscm_data(x2),'Vartype','unequal');
%     [tscmpu(1),tscmhu(1)] = ranksum(cd4tscm_data(x1),cd4tscm_data(x2));
% 
% exportgraphics(gcf,'response_CD4tscm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cd8tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd8tscm_data)*iqrmul,iqr(cd8tscm_data)*iqrmul];
% 
%     boxplot(cd8tscm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 Tscm')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd8tscm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd8tscm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8tscm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%    %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [tscmht(2),tscmpt(2)] = ttest2(cd8tscm_data(x1),cd8tscm_data(x2),'Vartype','unequal');
%     [tscmpu(2),tscmhu(2)] = ranksum(cd8tscm_data(x1),cd8tscm_data(x2));
% 
% exportgraphics(gcf,'response_CD8tscm.tiff','Resolution',1200)
% 
% % Responders / non-responders Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(tscm_data)*iqrmul,iqr(tscm_data)*iqrmul];
% 
%     boxplot(tscm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% Tscm')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(tscm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(tscm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,tscm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
% 
%     [tscmht(3),tscmpt(3)] = ttest2(tscm_data(x1),tscm_data(x2),'Vartype','unequal');
%     [tscmpu(3),tscmhu(3)] = ranksum(tscm_data(x1),tscm_data(x2));
% 
% 
% exportgraphics(gcf,'response_tscm.tiff','Resolution',1200)
%  
% % Parameter scatter Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
% m = m;
% n = n;
% t = tiledlayout(1,6);
% 
% for i = 1:size(parameters2,2)
% 
%     %subplot(m,n,i)
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     param = categorical(cellstr(par_names(i)));
% 
%     scatter(tscm_data, parameter_data,72,'kx','LineWidth',1.5);
% 
%     ylim(ylim_par(i,:));
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)  
%     xlabel('% Tscm')
%     ylabel(param)%,'Interpreter','none')
%     xlim([0 25])
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_Tscm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4 CM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd4cm_data)*iqrmul,iqr(cd4cm_data)*iqrmul];
% 
%     boxplot(cd4cm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 CM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd4cm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd4cm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4cm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
%     
%     [tcmht(1),tcmpt(1)] = ttest2(cd4cm_data(x1),cd4cm_data(x2),'Vartype','unequal');
%     [tcmpu(1),tcmhu(1)] = ranksum(cd4cm_data(x1),cd4cm_data(x2));
% 
% exportgraphics(gcf,'response_CD4cm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 Tcm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd8cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd8cm_data)*iqrmul,iqr(cd8cm_data)*iqrmul];
% 
%     boxplot(cd8cm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 CM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd8cm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd8cm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8cm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [tcmht(2),tcmpt(2)] = ttest2(cd8cm_data(x1),cd8cm_data(x2),'Vartype','unequal');
%     [tcmpu(2),tcmhu(2)] = ranksum(cd8cm_data(x1),cd8cm_data(x2));
% 
% exportgraphics(gcf,'response_CD8cm.tiff','Resolution',1200)
% 
% % Responders / non-responders TCM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cm_data)*iqrmul,iqr(cm_data)*iqrmul];
% 
%     boxplot(cm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [tcmht(3),tcmpt(3)] = ttest2(cm_data(x1),cm_data(x2),'Vartype','unequal');
%     [tcmpu(3),tcmhu(3)] = ranksum(cm_data(x1),cm_data(x2));
% 
% exportgraphics(gcf,'response_cm.tiff','Resolution',1200)
% 
% % Parameter scatter Tcm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
% m = m;
% n = n;
% t = tiledlayout(1,6);
% 
% 
% for i = 1:size(parameters2,2)
% 
%     %subplot(m,n,i)
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     param = categorical(cellstr(par_names(i)));
% 
%     scatter(cm_data, parameter_data,72,'kx','LineWidth',1.5);
% 
%     ylim(ylim_par(i,:));
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)  
%     xlabel('% CM')
%     xlim([0 100])
%     ylabel(param)%'Interpreter','none')
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_Tcm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4 EM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd4EM_data)*iqrmul,iqr(cd4EM_data)*iqrmul];
% 
%     boxplot(cd4EM_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 EM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd4EM_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd4EM_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4EM_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [temht(1),tempt(1)] = ttest2(cd4EM_data(x1),cd4EM_data(x2),'Vartype','unequal');
%     [tempu(1),temhu(1)] = ranksum(cd4EM_data(x1),cd4EM_data(x2));
% 
% exportgraphics(gcf,'response_CD4EM.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 TEM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd8EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd8EM_data)*iqrmul,iqr(cd8EM_data)*iqrmul];
% 
%     boxplot(cd8EM_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 EM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd8EM_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd8EM_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8EM_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [temht(2),tempt(2)] = ttest2(cd8EM_data(x1),cd8EM_data(x2),'Vartype','unequal');
%     [tempu(2),temhu(2)] = ranksum(cd8EM_data(x1),cd8EM_data(x2));
% 
% exportgraphics(gcf,'response_CD8EM.tiff','Resolution',1200)
% 
% % Responders / non-responders TEM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(EM_data)*iqrmul,iqr(EM_data)*iqrmul];
% 
%     boxplot(EM_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% EM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(EM_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(EM_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,EM_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [temht(3),tempt(3)] = ttest2(EM_data(x1),EM_data(x2),'Vartype','unequal');
%     [tempu(3),temhu(3)] = ranksum(EM_data(x1),EM_data(x2));
% 
% exportgraphics(gcf,'response_EM.tiff','Resolution',1200)
% 
% 
% % Parameter scatter TsEM ratio
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width,par_plot_height])
% m = m;
% n = n;
% t = tiledlayout(1,6);
% 
% for i = 1:size(parameters2,2)
% 
%     %subplot(m,n,i)
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     param = categorical(cellstr(par_names(i)));
% 
%     scatter(EM_data, parameter_data,72,'kx','LineWidth',1.5);
% 
%     ylim(ylim_par(i,:));
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)  
%     xlabel('% EM')
%     xlim([0 100])
%     ylabel(param)%,'Interpreter','none')
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_TEM.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4 EMRA
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd4EMRA_data)*iqrmul,iqr(cd4EMRA_data)*iqrmul];
% 
%     boxplot(cd4EMRA_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 EMRA')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd4EMRA_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd4EMRA_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4EMRA_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [temraht(1),temrapt(1)] = ttest2(cd4EMRA_data(x1),cd4EMRA_data(x2),'Vartype','unequal');
%     [temrapu(1),temrahu(1)] = ranksum(cd4EMRA_data(x1),cd4EMRA_data(x2));
% 
% exportgraphics(gcf,'response_CD4EMRA.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 TEMRA
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd8EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(cd8EMRA_data)*iqrmul,iqr(cd8EMRA_data)*iqrmul];
% 
%     boxplot(cd8EMRA_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 EMRA')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(cd8EMRA_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(cd8EMRA_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8EMRA_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [temraht(2),temrapt(2)] = ttest2(cd8EMRA_data(x1),cd8EMRA_data(x2),'Vartype','unequal');
%     [temrapu(2),temrahu(2)] = ranksum(cd8EMRA_data(x1),cd8EMRA_data(x2));
% 
% exportgraphics(gcf,'response_CD8EMRA.tiff','Resolution',1200)
% 
% % Responders / non-responders TEMRA
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Day28_response2')));
%     DL = [-iqr(EMRA_data)*iqrmul,iqr(EMRA_data)*iqrmul];
% 
%     boxplot(EMRA_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'Response', 'No response'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% EMRA')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'Response'); x(x1) = 1; par_mean_x1 = mean(EMRA_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'No response'); x(x2) = 2; par_mean_x2 = mean(EMRA_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,EMRA_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'R','NR'})
% 
%     [temraht(3),temrapt(3)] = ttest2(EMRA_data(x1),EMRA_data(x2),'Vartype','unequal');
%     [temrapu(3),temrahu(3)] = ranksum(EMRA_data(x1),EMRA_data(x2));
% 
% exportgraphics(gcf,'response_EMRA.tiff','Resolution',1200)
% 
% % Parameter scatter TsEMRA ratio
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width,par_plot_height])
% m = m;
% n = n;
% t = tiledlayout(1,6);
% 
% 
% for i = 1:size(parameters2,2)
% 
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     param = categorical(cellstr(par_names(i)));
% 
%     scatter(EMRA_data, parameter_data,72,'kx','LineWidth',1.5);
% 
%     ylim(ylim_par(i,:));
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)  
%     xlabel('% EMRA')
%     xlim([0 100])
%     ylabel(param)%,'Interpreter','none')
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_TEMRA.tiff','Resolution',1200)
% 
% 
% % Parameter scatter TE%
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width,par_plot_height])
% m = m;
% n = n;
% t = tiledlayout(1,6);
% 
% for i = 1:size(parameters2,2)
% 
%     %subplot(m,n,i)
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     TE_data = (table2array(var(:,'% CAR+ of CD3+ (aka transduction efficiency):')))*100;
%     param = categorical(cellstr(par_names(i)));
% 
%     scatter(TE_data, parameter_data,72,'kx','LineWidth',1.5);
% 
%     ylim(ylim_par(i,:));
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)  
%     xlabel('% TE')
%     xlim([0 100])
%     ylabel(param)%,'Interpreter','none')
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_TE.tiff','Resolution',1200)


% %% Parameter scatter disease
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/2,par_plot_height*2])
% %sgtitle("Disease")
% t = tiledlayout(2,3);
% 
% for i = 1:size(parameters2,2)
% 
%     %subplot(m,n,i)
%     nexttile
%     parameter_data = table2array(var(:,parameters2(i)));
%     categorical_data = categorical(table2array(var(:,'Disease')));
%     param = categorical(cellstr(par_names(i)));
%     DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];
% 
%     boxplot(parameter_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'DLBCL','MCL', 'FL', 'CLL'}); 
% 
%     ax = gca;
%     ax.YAxis.Exponent = y_exp_par(i);
%     set(gca,'XTickLabel',[]);
%     set(gca, 'YScale', par_yscale(i))
%     set(gca,'FontSize',par_fontsize)
%     ylim(ylim_par(i,:));
%     xlabel(param)%,'Interpreter','none')
%     xtickangle(45)
%     
%     hold on;
%     x = zeros([size(categorical_data,1),1]);
%     x1 = (categorical_data == 'DLBCL'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'MCL'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
%     x3 = (categorical_data == 'FL'); x(x3) = 3; par_mean_x3 = mean(parameter_data(x3)); par_n_x3 = sum(x3);
%     x4 = (categorical_data == 'CLL'); x(x4) = 4; par_mean_x4 = mean(parameter_data(x4)); par_n_x4 = sum(x4);
% 
%     swarmchart(x(x1),parameter_data(x1),[],'g');
%     swarmchart(x(x2),parameter_data(x2),[],'MarkerEdgeColor','#EDB120');
%     swarmchart(x(x3),parameter_data(x3),[],'c');
%     swarmchart(x(x4),parameter_data(x4),[],'m');
%     plot(1,par_mean_x1,'og','MarkerSize',7,'MarkerFaceColor','g');
%     plot(2,par_mean_x2,'o','MarkerSize',7,'MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120');
%     plot(3,par_mean_x3,'oc','MarkerSize',7,'MarkerFaceColor','c');
%     plot(4,par_mean_x4,'om','MarkerSize',7,'MarkerFaceColor','m');
% 
%     if i == 8
%         xticklabels({['DLBCL (n=',num2str(par_n_x1),')'],['MCL (n=',num2str(par_n_x2),')'], ...
%                     ['FL (n=',num2str(par_n_x3),')'],['CLL (n=',num2str(par_n_x4),')']})
%     else
%         xticklabels({'DLBCL','MCL','FL','CLL'})
%     end
% 
% %     x = zeros(size(categorical_data,1));
% %     x1 = (categorical_data == 'DLBCL'); x(x1) = 1;
% %     x2 = (categorical_data == 'MCL'); x(x2) = 2;
% %     x3 = (categorical_data == 'FL'); x(x3) = 3;
% %     x4 = (categorical_data == 'CLL'); x(x4) = 4;
% %     
% %     swarmchart(x,parameter_data,'xb');
% 
%     %p_disease(i) = kruskalwallis(parameter_data, categorical_data);
% 
%     hold off;
% 
% end
% 
% exportgraphics(gcf,'parameters_dist_disease.tiff','Resolution',1200)


%% Relapse

included_patients = var(strcmp(var.Day28_response2, 'Response'),'Patient');
[patients,ilp,irp] = innerjoin(patients,included_patients,"LeftKeys","patients","RightKeys","Patient");
time = time(:,ilp);
tumor = tumor(:,ilp);
e25t1 = e25t1(:,ilp);
e6t1 = e6t1(:,ilp);
e1t1 = e1t1(:,ilp);
cytolysis = cytolysis(:,ilp);
e25 = e25(:,ilp);
e6 = e6(:,ilp);
e1 = e1(:,ilp);
Ein = Ein(:,ilp);
e25t1_sim = e25t1_sim(:,ilp);
e6t1_sim = e6t1_sim(:,ilp);
e1t1_sim = e1t1_sim(:,ilp);
tumor_sim = tumor_sim(:,ilp);
var = var(strcmp(var.Day28_response2, 'Response'),:);

% Parameter scatter relapse
m = 1;
n = 9;

f = f + 1;
figure(f)
set(gcf,'Units','inches','Position',[0.5,0.5,par_plot_width/1.5,par_plot_height])
%sgtitle("Relapse")
t = tiledlayout(1,6);

for i = 1:size(parameters2,2)

    %subplot(m,n,i)
    nexttile
    parameter_data = table2array(var(:,parameters2(i)));
    categorical_data = categorical(table2array(var(:,'Relapse')));
    param = categorical(cellstr(par_names(i)));
    DL = [-iqr(parameter_data)*iqrmul,iqr(parameter_data)*iqrmul];

    boxplot(parameter_data, categorical_data,...
            DataLim=DL,...
            ExtremeMode="clip",...
            Colors='k',...
            PlotStyle='traditional',...
            OutlierSize=0.1,...
            Widths=0.5,...
            GroupOrder={'No','Yes'}); 

    ax = gca;
    ax.YAxis.Exponent = y_exp_par(i);
    set(gca,'XTickLabel',[]);
    set(gca, 'YScale', par_yscale(i))
    set(gca,'FontSize',par_fontsize)
    ylim(ylim_par(i,:));
    xlabel(param)%,'Interpreter','none')
    
    hold on;
    x = zeros([size(categorical_data,1),1]);
    x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(parameter_data(x1)); par_n_x1 = sum(x1);
    x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(parameter_data(x2)); par_n_x2 = sum(x2);
    
    swarmchart(x(x1),parameter_data(x1),[],[0 0 1]);
    swarmchart(x(x2),parameter_data(x2),[],[1 0 0]);
    plot(1,par_mean_x1,'ob','MarkerSize',7,'MarkerFaceColor','b');
    plot(2,par_mean_x2,'or','MarkerSize',7,'MarkerFaceColor','r');

    if i == 8
        xticklabels({['No (n=',num2str(par_n_x2),')'],['Yes (n=',num2str(par_n_x1),')']})
    else
        xticklabels({'N','Y'})
    end

    % statistics for each group
    [chu(i,3),cpu(i,3),stats] = ranksum(parameter_data(x1),parameter_data(x2));
    w = stats.ranksum;
    q = (sum(x2)*(sum(x2)+1))/2;
    u = w - q;
    r = 1 - ((2*u)/(sum(x1)*sum(x2)));
    ci_relapse(i,1) = r;
    medDiff = @(x,y) median(x) - median(y);
    stats_relapse_CI(i,:) = meanEffectSize(parameter_data(x1),parameter_data(x2),...
                        "VarianceType","unequal",ConfidenceIntervalType="bootstrap");


    hold off;

end

exportgraphics(gcf,'parameters_dist_relapse.tiff','Resolution',1200)

% % Relapse
% 
% % data CR
% patients_CR = var(strcmp(var.Relapse,'No'),"Patient");
% [patients_CR,ilpcr,irpcr] = innerjoin(patients,patients_CR,"LeftKeys","patients","RightKeys","Patient");
% time_CR = time(:,ilpcr);
% tumor_CR = tumor(:,ilpcr);
% e25t1_CR = e25t1(:,ilpcr);
% e6t1_CR = e6t1(:,ilpcr);
% e1t1_CR = e1t1(:,ilpcr);
% 
% % means CR
% [std_time_CR,mean_time_CR] = std(time_CR,0,2,"omitnan");
% [std_tumor_CR,mean_tumor_CR] = std(tumor_CR,0,2,"omitnan");
% [std_e25t1_CR,mean_e25t1_CR] = std(e25t1_CR,0,2,"omitnan");
% [std_e6t1_CR,mean_e6t1_CR] = std(e6t1_CR,0,2,"omitnan");
% [std_e1t1_CR,mean_e1t1_CR] = std(e1t1_CR,0,2,"omitnan");
% 
% % data NCR
% patients_NCR = var(strcmp(var.Relapse,'Yes'),"Patient");
% [patients_NCR,ilpncr,irpncr] = innerjoin(patients,patients_NCR,"LeftKeys","patients","RightKeys","Patient");
% time_NCR = time(:,ilpncr);
% tumor_NCR = tumor(:,ilpncr);
% e25t1_NCR = e25t1(:,ilpncr);
% e6t1_NCR = e6t1(:,ilpncr);
% e1t1_NCR = e1t1(:,ilpncr);
% 
% % means NCR
% [std_time_NCR,mean_time_NCR] = std(time_NCR,0,2,"omitnan");
% [std_tumor_NCR,mean_tumor_NCR] = std(tumor_NCR,0,2,"omitnan");
% [std_e25t1_NCR,mean_e25t1_NCR] = std(e25t1_NCR,0,2,"omitnan");
% [std_e6t1_NCR,mean_e6t1_NCR] = std(e6t1_NCR,0,2,"omitnan");
% [std_e1t1_NCR,mean_e1t1_NCR] = std(e1t1_NCR,0,2,"omitnan");
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,16,4])
% startColor = [0.5, 0.5, 0.5];
% %sgtitle("cytotoxicity by day 28 response (blue - CR, red - NCR)")
% 
% subplot(1,4,1)
% % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e1t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e1t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 1:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,2)
% % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e6t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e6t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 6.25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(1,4,3)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     x1 = (categorical_data == 'No'); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); par_n_x2 = sum(x2);
% 
% subplot(1,4,4)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% axis([10,11,10,11])
% H = legend( ['No Relapse (N) (n=',num2str(par_n_x1),')'],...
%             'SD N',...
%             ['Relapse (Y) (n=',num2str(par_n_x2),')'],...
%             'SD Y',...
%             'Location','west');
% set(H,'FontSize',legend_fontsize)
% legend boxoff
% axis off
% 
% % each time course
% 
% % subplot(3,3,4)
% % 
% % plot(mean_time_CR,mean_e25t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e25t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,7)
% % 
% % plot(mean_time_NCR,mean_e25t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e25t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,5)
% % 
% % plot(mean_time_CR,mean_e6t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e6t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,8)
% % 
% % plot(mean_time_NCR,mean_e6t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e6t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,6)
% % 
% % plot(mean_time_CR,mean_e1t1_CR,'--b','Linewidth',2.5); hold on;
% % patch([mean_time_CR' flip(mean_time_CR')],...
% %       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
% %       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_CR,2)
% %     plot(time_CR(:,ii),e1t1_CR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% % subplot(3,3,9)
% % 
% % plot(mean_time_NCR,mean_e1t1_NCR,'--r','Linewidth',2.5); hold on;
% % patch([mean_time_NCR' flip(mean_time_NCR')],...
% %       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
% %       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% % for ii = 1:size(time_NCR,2)
% %     plot(time_NCR(:,ii),e1t1_NCR(:,ii),'Color', startColor + ii/100);
% % end
% % xlim([0 72]);
% % xlabel('Time (hours)'); ylabel('CI');
% % set(gca,'FontSize',16)
% % 
% 
% exportgraphics(gcf,'all_data_relapse.tiff','Resolution',1200)
% 
% % Relapse each time course
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,12,12])
% startColor = [0.5, 0.5, 0.5];
% %sgtitle("cytotoxicity by relapse (blue - CR, red - NCR)")
% m = 3;
% n = 3;
% 
% subplot(m,n,1)
% % plot(mean_time,mean_e1t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e1t1 + std_e1t1)' flip((mean_e1t1 - std_e1t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e1t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e1t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 1:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(m,n,2)
% % plot(mean_time,mean_e6t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e6t1 + std_e6t1)' flip((mean_e6t1 - std_e6t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e6t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e6t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 6.25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% subplot(m,n,3)
% % plot(mean_time,mean_e25t1,'--k','Linewidth',2.5); hold on;
% % patch([mean_time' flip(mean_time')],...
% %       [(mean_e25t1 + std_e25t1)' flip((mean_e25t1 - std_e25t1)')],...
% %       [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_CR,mean_e25t1_CR,'-b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% plot(mean_time_NCR,mean_e25t1_NCR,'-r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% title("E:T 25:1")
% % H = legend('all','SD all', 'Response',...
% %             'SD Response','No response', 'SD No response',...
% %             'Location','northeastoutside');
% % set(H,'FontSize',legend_fontsize)
% % legend boxoff
% 
% %each time course
% 
% startColor = [0.0, 0.0, 0.5];
% 
% subplot(m,n,6)
% 
% plot(mean_time_CR,mean_e25t1_CR,'--b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e25t1_CR + std_e25t1_CR)' flip((mean_e25t1_CR - std_e25t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_CR,2)
%     plot(time_CR(:,ii),e25t1_CR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,5)
% 
% plot(mean_time_CR,mean_e6t1_CR,'--b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e6t1_CR + std_e6t1_CR)' flip((mean_e6t1_CR - std_e6t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_CR,2)
%     plot(time_CR(:,ii),e6t1_CR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,4)
% 
% plot(mean_time_CR,mean_e1t1_CR,'--b','Linewidth',2.5); hold on;
% patch([mean_time_CR' flip(mean_time_CR')],...
%       [(mean_e1t1_CR + std_e1t1_CR)' flip((mean_e1t1_CR - std_e1t1_CR)')],...
%       [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_CR,2)
%     plot(time_CR(:,ii),e1t1_CR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% startColor = [0.5, 0.0, 0.0];
% 
% subplot(m,n,9)
% 
% plot(mean_time_NCR,mean_e25t1_NCR,'--r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e25t1_NCR + std_e25t1_NCR)' flip((mean_e25t1_NCR - std_e25t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_NCR,2)
%     plot(time_NCR(:,ii),e25t1_NCR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,8)
% 
% plot(mean_time_NCR,mean_e6t1_NCR,'--r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e6t1_NCR + std_e6t1_NCR)' flip((mean_e6t1_NCR - std_e6t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_NCR,2)
%     plot(time_NCR(:,ii),e6t1_NCR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% subplot(m,n,7)
% 
% plot(mean_time_NCR,mean_e1t1_NCR,'--r','Linewidth',2.5); hold on;
% patch([mean_time_NCR' flip(mean_time_NCR')],...
%       [(mean_e1t1_NCR + std_e1t1_NCR)' flip((mean_e1t1_NCR - std_e1t1_NCR)')],...
%       [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
% for ii = 1:size(time_NCR,2)
%     plot(time_NCR(:,ii),e1t1_NCR(:,ii),'Color', startColor + ii/100);
% end
% xlim(xlim_time);
% ylim(ylim_CI);
% xticks(xticks_time);
% yticks(yticks_CI);
% xlabel('Time (hours)'); ylabel('CI');
% set(gca,'FontSize',plot_fontsize2)
% 
% exportgraphics(gcf,'all_data_relapse_tc.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4:CD8 ratio
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cd4cd8_data = table2array(var(:,'CD4_CD8_Ratio_CART'));
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd4cd8_data)*iqrmul,iqr(cd4cd8_data)*iqrmul];
% 
%     boxplot(cd4cd8_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('CD4:CD8')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd4cd8_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd4cd8_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4cd8_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     xticklabels({['N (n=',num2str(par_n_x1),')'],['Y (n=',num2str(par_n_x2),')']})
%     xtickangle(0)
% 
%     % statistics for each group
%     [cd4cd8ht(3),cd4cd8pt(3)] = ttest2(cd4cd8_data(x1),cd4cd8_data(x2),'Vartype','unequal');
%     [cd4cd8pu(3),cd4cd8hu(3)] = ranksum(cd4cd8_data(x1),cd4cd8_data(x2));
% 
% exportgraphics(gcf,'relapse_CD4CD8.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4 Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd4tscm_data)*iqrmul,iqr(cd4tscm_data)*iqrmul];
% 
%     boxplot(cd4tscm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 Tscm')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd4tscm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd4tscm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4tscm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'N','Y'})
% 
%     [tscmht_re(1),tscmpt_re(1)] = ttest2(cd4tscm_data(x1),cd4tscm_data(x2),'Vartype','unequal');
%     [tscmpu_re(1),tscmhu_re(1)] = ranksum(cd4tscm_data(x1),cd4tscm_data(x2));
% 
% exportgraphics(gcf,'relapse_CD4tscm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cd8tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd8tscm_data)*iqrmul,iqr(cd8tscm_data)*iqrmul];
% 
%     boxplot(cd8tscm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 Tscm')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd8tscm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd8tscm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8tscm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%    %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [tscmht_re(2),tscmpt_re(2)] = ttest2(cd8tscm_data(x1),cd8tscm_data(x2),'Vartype','unequal');
%     [tscmpu_re(2),tscmhu_re(2)] = ranksum(cd8tscm_data(x1),cd8tscm_data(x2));
% 
% exportgraphics(gcf,'relapse_CD8tscm.tiff','Resolution',1200)
% 
% % Responders / non-responders Tscm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     tscm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8Tscm')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(tscm_data)*iqrmul,iqr(tscm_data)*iqrmul];
% 
%     boxplot(tscm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% Tscm')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(tscm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(tscm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,tscm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     xticklabels({['No (n=',num2str(par_n_x1),')'],['Yes (n=',num2str(par_n_x2),')']})
% 
%     [tscmht_re(3),tscmpt_re(3)] = ttest2(tscm_data(x1),tscm_data(x2),'Vartype','unequal');
%     [tscmpu_re(3),tscmhu_re(3)] = ranksum(tscm_data(x1),tscm_data(x2));
% 
% 
% exportgraphics(gcf,'relapse_tscm.tiff','Resolution',1200)
%  
% % Responders / non-responders CD4 CM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd4cm_data)*iqrmul,iqr(cd4cm_data)*iqrmul];
% 
%     boxplot(cd4cm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 CM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd4cm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd4cm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4cm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
%     
%     [tcmht_re(1),tcmpt_re(1)] = ttest2(cd4cm_data(x1),cd4cm_data(x2),'Vartype','unequal');
%     [tcmpu_re(1),tcmhu_re(1)] = ranksum(cd4cm_data(x1),cd4cm_data(x2));
% 
% exportgraphics(gcf,'relapse_CD4cm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 Tcm
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd8cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd8cm_data)*iqrmul,iqr(cd8cm_data)*iqrmul];
% 
%     boxplot(cd8cm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 CM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd8cm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd8cm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8cm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [tcmht_re(2),tcmpt_re(2)] = ttest2(cd8cm_data(x1),cd8cm_data(x2),'Vartype','unequal');
%     [tcmpu_re(2),tcmhu_re(2)] = ranksum(cd8cm_data(x1),cd8cm_data(x2));
% 
% exportgraphics(gcf,'relapse_CD8cm.tiff','Resolution',1200)
% 
% % Responders / non-responders TCM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
% %sgtitle("CD4:CD8 ratio")
%     
%     cm_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8CM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cm_data)*iqrmul,iqr(cm_data)*iqrmul];
% 
%     boxplot(cm_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cm_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cm_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cm_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [tcmht_re(3),tcmpt_re(3)] = ttest2(cm_data(x1),cm_data(x2),'Vartype','unequal');
%     [tcmpu_re(3),tcmhu_re(3)] = ranksum(cm_data(x1),cm_data(x2));
% 
% exportgraphics(gcf,'relapse_cm.tiff','Resolution',1200)
% 
% % Responders / non-responders CD4 EM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd4EM_data)*iqrmul,iqr(cd4EM_data)*iqrmul];
% 
%     boxplot(cd4EM_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 EM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd4EM_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd4EM_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4EM_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [temht_re(1),tempt_re(1)] = ttest2(cd4EM_data(x1),cd4EM_data(x2),'Vartype','unequal');
%     [tempu_re(1),temhu_re(1)] = ranksum(cd4EM_data(x1),cd4EM_data(x2));
% 
% exportgraphics(gcf,'relapse_CD4EM.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 TEM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd8EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd8EM_data)*iqrmul,iqr(cd8EM_data)*iqrmul];
% 
%     boxplot(cd8EM_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 EM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd8EM_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd8EM_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8EM_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [temht_re(2),tempt_re(2)] = ttest2(cd8EM_data(x1),cd8EM_data(x2),'Vartype','unequal');
%     [tempu_re(2),temhu_re(2)] = ranksum(cd8EM_data(x1),cd8EM_data(x2));
% 
% exportgraphics(gcf,'relapse_CD8EM.tiff','Resolution',1200)
% 
% % Responders / non-responders TEM
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     EM_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EM')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(EM_data)*iqrmul,iqr(EM_data)*iqrmul];
% 
%     boxplot(EM_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% EM')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(EM_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(EM_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,EM_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [temht_re(3),tempt_re(3)] = ttest2(EM_data(x1),EM_data(x2),'Vartype','unequal');
%     [tempu_re(3),temhu_re(3)] = ranksum(EM_data(x1),EM_data(x2));
% 
% exportgraphics(gcf,'relapse_EM.tiff','Resolution',1200)
% 
% 
% % Responders / non-responders CD4 EMRA
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd4EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd4EMRA_data)*iqrmul,iqr(cd4EMRA_data)*iqrmul];
% 
%     boxplot(cd4EMRA_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD4 EMRA')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd4EMRA_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd4EMRA_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd4EMRA_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [temraht_re(1),temrapt_re(1)] = ttest2(cd4EMRA_data(x1),cd4EMRA_data(x2),'Vartype','unequal');
%     [temrapu_re(1),temrahu_re(1)] = ranksum(cd4EMRA_data(x1),cd4EMRA_data(x2));
% 
% exportgraphics(gcf,'relapse_CD4EMRA.tiff','Resolution',1200)
% 
% % Responders / non-responders CD8 TEMRA
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     cd8EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(cd8EMRA_data)*iqrmul,iqr(cd8EMRA_data)*iqrmul];
% 
%     boxplot(cd8EMRA_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% CD8 EMRA')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(cd8EMRA_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(cd8EMRA_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,cd8EMRA_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [temraht_re(2),temrapt_re(2)] = ttest2(cd8EMRA_data(x1),cd8EMRA_data(x2),'Vartype','unequal');
%     [temrapu_re(2),temrahu_re(2)] = ranksum(cd8EMRA_data(x1),cd8EMRA_data(x2));
% 
% exportgraphics(gcf,'relapse_CD8EMRA.tiff','Resolution',1200)
% 
% % Responders / non-responders TEMRA
% 
% f = f + 1;
% figure(f)
% set(gcf,'Units','inches','Position',[0.5,0.5,3,par_plot_height])
%     
%     EMRA_data = (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD4EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD4')))*100 + ...
%                     (table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD8EMRA')) .* ...
%                     table2array(var(:,'NEWDIFFERENTIATIONPANEL_CD3CD8')))*100;
%     categorical_data = categorical(table2array(var(:,'Relapse')));
%     DL = [-iqr(EMRA_data)*iqrmul,iqr(EMRA_data)*iqrmul];
% 
%     boxplot(EMRA_data, categorical_data,...
%             DataLim=DL,...
%             ExtremeMode="clip",...
%             Colors='k',...
%             PlotStyle='traditional',...
%             OutlierSize=0.1,...
%             Widths=0.5,...
%             GroupOrder={'No', 'Yes'}); 
% 
%     %ylim([0 20]);
%     set(gca,'XTickLabel',[]);
%     set(gca,'FontSize',par_fontsize)
%     ylabel('% EMRA')
%     
%     hold on;
%     x = zeros(size(categorical_data,1));
%     x1 = (categorical_data == 'No'); x(x1) = 1; par_mean_x1 = mean(EMRA_data(x1)); par_n_x1 = sum(x1);
%     x2 = (categorical_data == 'Yes'); x(x2) = 2; par_mean_x2 = mean(EMRA_data(x2)); par_n_x2 = sum(x2);
%     
%     swarmchart(x,EMRA_data,[],[0.5 0.5 0.5]);
%     plot(1,par_mean_x1,'ok','MarkerSize',7,'MarkerFaceColor','k');
%     plot(2,par_mean_x2,'ok','MarkerSize',7,'MarkerFaceColor','k');
% 
%     %xticklabels({['R (n=',num2str(par_n_x1),')'],['NR (n=',num2str(par_n_x2),')']})
%     xticklabels({'No','Yes'})
% 
%     [temraht_re(3),temrapt_re(3)] = ttest2(EMRA_data(x1),EMRA_data(x2),'Vartype','unequal');
%     [temrapu_re(3),temrahu_re(3)] = ranksum(EMRA_data(x1),EMRA_data(x2));
% 
% exportgraphics(gcf,'relapse_EMRA.tiff','Resolution',1200)
