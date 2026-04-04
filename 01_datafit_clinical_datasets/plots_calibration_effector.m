function f = plots_calibration_effector(CI, cells, qfit)

cells_range = (0:1000:2000000);
CI_calc = polyval(qfit,cells_range);

f = figure("Visible","off");

plot(CI,cells,'*r','MarkerSize',16); hold on;
plot(CI_calc, cells_range, '--k','LineWidth',3);
ylim([0 2000000]);
xlim([-0.1 2])
xlabel('CI'); ylabel('Cells');
set(gca,'Fontsize',20); 
H = legend('Effector Data',' Convesion (cells->CI)', 'Location','southeast');
set(H,'FontSize',16)
legend boxoff

end