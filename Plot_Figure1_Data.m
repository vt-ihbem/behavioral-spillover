%% Load and Plot Data of COVID and Flu cases
clear; close all; clc

saveplot = 1;

%% Load Data 

COVID_data = readtable('Fig 1 COVID-Flu-US Data.xlsx','Sheet','COVID-19 cases');
Flu_data = readtable('Fig 1 COVID-Flu-US Data.xlsx','Sheet','Flu');

COVID_names = COVID_data.Properties.VariableNames;
Flu_names = Flu_data.Properties.VariableNames;

%% Plot Data
f1 = figure(1);
clf
yyaxis left
plot(COVID_data.Day,COVID_data.WeeklyCOVID_19Cases,'LineWidth',2)
ylabel('COVID-19 cases (thousands)')
set(gca,'yticklabel',0:1000:6000)
hold on
yyaxis right
plot(Flu_data.Day,Flu_data.InfluenzaA_AllTypesOfSurveillance,'LineWidth',2)
plot(Flu_data.Day,Flu_data.InfluenzaB_AllTypesOfSurveillance,'LineWidth',2)
hold off
set(gca,'fontsize',18)
ylabel('Influenza cases (thousands)')
set(gca,'yticklabel',0:10:60)
l = legend('COVID-19','Influenza A','Influenza B');
set(l,'box','off','location','northwest')
xlim([datetime(2016,01,01) datetime(2025,01,01)])
xtickangle(45)

%% Save plot

if saveplot == 1
    figname = 'Figures/COVIDFluData.png';
    saveas(f1,figname,'png')
    figname1 = 'Figures/tiff/COVIDFluData..tif';
    saveas(f1,figname1,'tif')
end
