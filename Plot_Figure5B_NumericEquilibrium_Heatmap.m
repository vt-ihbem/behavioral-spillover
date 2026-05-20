%% Heatmap of threshold
clearvars; close all; clc

saveplot = 1; % 1 to save plots

filein = 'Data_Numeric_Equilibrium.mat';
load(filein)

%%
f1 = figure(1);
h = imagesc(sin_vec,R02_vec,thres>R0out);
set(h,'AlphaData',0.75)
C=jet(4);
C(1,:) = [0.4940, 0.1840, 0.5560];
colormap(C)
set(gca, 'Ydir','normal','fontsize',24,'fontname','Times')
xlabel('spillover, $s$','Interpreter','latex')
ylabel('$\mathcal{R}_{0,B}$','Interpreter','latex')
axis square
hold on
s = [0 .5 1 0 .5 1 0 .5 1];
R = [2.9 2.9 2.9 2 2 2 1.3 1.3 1.3];
plot(s,R,'.','color',[.8 .8 .8],'markersize',40,'linewidth',2)
plot(s,R,'k.','markersize',20)
text(.2,2.5,'Co-existence','fontsize',22,'interpreter','latex','color','w')
text(.58,1.22,'Exclusion','fontsize',22,'interpreter','latex','color','w')
text(.56,1.1,'of disease $B$','fontsize',22,'interpreter','latex','color','w')

%% Save Figure
if saveplot==1
    % png
    fileout = strcat('Figures/Heatmap_NumericEquilibrium.png');
    saveas(f1,fileout,'png')
    % tif
    fileout1 = strcat('Figures/tiff/Heatmap_NumericEquilibrium.tif');
    saveas(f1,fileout1,'tif')
end
