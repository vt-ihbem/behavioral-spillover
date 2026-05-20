%% Plot of Approximation Threshold
clearvars; clc; close all; clc

saveplot = 1;
R0B_vec = 1:.02:3;

R0A = 3;

sthres = NaN(length(R0B_vec),1);
for it = 1:length(R0B_vec)
    R0B = R0B_vec(it);

    sthres(it) = (1-1/R0B)/(1-1/R0A);
end

f1 = figure(1);
plot(sthres,R0B_vec,'r','linewidth',2)
xlabel('spillover, $s$','Interpreter','latex')
ylabel('$\mathcal{R}_{0,B}$','Interpreter','latex')
set(gca,'fontsize',24,'fontname','times')
text(.5,1.58,'$s_{threshold}$','fontsize',16,'Interpreter','latex','Rotation',45,'color','r')
text(.2,2.5,'Co-existence','fontsize',20,'interpreter','latex')
text(.45,1.3,'Potential exclusion','fontsize',20,'interpreter','latex')
text(.55,1.18,'of disease $B$','fontsize',20,'interpreter','latex')
axis square

%% Save Figure
if saveplot == 1
    % png
    fileout = strcat('Figures/Threshold.png');
    saveas(f1,fileout,'png')
    % tif
    fileouta = strcat('Figures/tiff/Threshold.tif');
    saveas(f1,fileouta,'tif')
end
