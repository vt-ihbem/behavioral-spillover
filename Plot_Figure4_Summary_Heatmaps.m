%% Load Data
clearvars; close all; clc

saveplot = 1; % 1 to save plots
files = {'Data_mu0_tauR100_long'};

figure(5); clf
figure(6); clf

filein = strcat(files{1},'.mat');
load(filein)

%% Summary Heatmap Plots

% Numerical simulation for co-existence vs exclusion
f5 = figure(5);
thres = 1e-5;
h = imagesc(svec,R02vec,iout2<thres);
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
hold off


% Time when disease B > disease A
Bdom3 = Bdom;
Bdom3(end,1) = NaN;
Bdom3(iout2<thres) = NaN;

f6 = figure(6);
h = imagesc(svec,R02vec,Bdom3*100);
set(h,'AlphaData', ~isnan(Bdom3))
set(gca, 'Ydir','normal','fontsize',24)
set(gca,'fontname','Times')
xlabel('spillover, $s$','Interpreter','latex')
ylabel('$\mathcal{R}_{0,B}$','Interpreter','latex')
c = colorbar;
c.Label.String = {'Percentage of time B exceeds A ',sprintf('(during t=0 to t=%d)',tmax)};
c.Label.FontName = 'Times';
axis square
hold on
s = [0 .5 1 0 .5 1 0 .5 1];
R = [2.9 2.9 2.9 2 2 2 1.3 1.3 1.3];
plot(s,R,'.','color',[.8 .8 .8],'markersize',40,'linewidth',2)
plot(s,R,'k.','markersize',20)
hold off


% Amount that disease B exceeds A
Bdom3 = Bdom1;
Bdom3(end,1) = NaN;
Bdom3(iout2<thres) = NaN;

f7 = figure(7);
h = imagesc(svec,R02vec,(Bdom3/3.7187*100));
set(h,'AlphaData', ~isnan(Bdom3))
set(gca, 'Ydir','normal','fontsize',24)
set(gca,'fontname','Times')
xlabel('spillover, $s$','Interpreter','latex')
ylabel('$\mathcal{R}_{0,B}$','Interpreter','latex')
c = colorbar;
c.Label.String = {'Cumulative amount by which B exceeds ',sprintf('(when it exceeds A in t=0 to t=%d)',tmax)};
c.Label.FontName = 'Times';
axis square
hold on
s = [0 .5 1 0 .5 1 0 .5 1];
R = [2.9 2.9 2.9 2 2 2 1.3 1.3 1.3];
plot(s,R,'.','color',[.8 .8 .8],'markersize',40,'linewidth',2)
plot(s,R,'k.','markersize',20)
hold off

% Amount that disease A exceeds B
Bdom3 = Bdom1a;
Bdom3(end,1) = NaN;
Bdom3(iout2<thres) = NaN;

f8 = figure(8);
h = imagesc(svec,R02vec,(Bdom3/3.7187*100));
set(h,'AlphaData', ~isnan(Bdom3))
set(gca, 'Ydir','normal','fontsize',24)
set(gca,'fontname','Times')
xlabel('spillover, $s$','Interpreter','latex')
ylabel('$\mathcal{R}_{0,B}$','Interpreter','latex')
c = colorbar;
c.Label.String = {'Cumulative amount by which A exceeds ',sprintf('(when it exceeds B in t=0 to t=%d)',tmax)};
c.Label.FontName = 'Times';
axis square
hold on
s = [0 .5 1 0 .5 1 0 .5 1];
R = [2.9 2.9 2.9 2 2 2 1.3 1.3 1.3];
plot(s,R,'.','color',[.8 .8 .8],'markersize',40,'linewidth',2)
plot(s,R,'k.','markersize',20)
hold off

%% Save Figures
if saveplot==1
    % png figures
    fileout5 = strcat('Figures/Heatmap_','Ithres.png');
    saveas(f5,fileout5,'png')
    fileout6 = strcat('Figures/Heatmap_','BoverA_time.png');
    saveas(f6,fileout6,'png')
    fileout7 = strcat('Figures/Heatmap_','BoverA.png');
    saveas(f7,fileout7,'png')
    fileout8 = strcat('Figures/Heatmap_','AoverB.png');
    saveas(f8,fileout8,'png')
    % tif figures
    fileout5a = strcat('Figures/tiff/Heatmap_','Ithres.tif');
    saveas(f5,fileout5a,'tif')
    fileout6a = strcat('Figures/tiff/Heatmap_','BoverA_time.tif');
    saveas(f6,fileout6a,'tif')
    fileout7a = strcat('Figures/tiff/Heatmap_','BoverA.tif');
    saveas(f7,fileout7a,'tif')
    fileout8a = strcat('Figures/tiff/Heatmap_','AoverB.tif');
    saveas(f8,fileout8a,'tif')
end


