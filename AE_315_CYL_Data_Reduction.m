% Code Name: CYL Aerodynamics Lab Data Reductor
% Code Description: Reduces the data collected in the CYL lab and plots
% necessary graphs and tables
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE315 - Section 29DB
% Date: 02/17/2026
% Worked With: N/a

% Iniciation
clear; clc; close all;
data15_t = load("29DB_Cyl_15mpers2_20260210T173458.mat");
data15 = data15_t.rawData;
data10_t = load('29DB_Cyl_10mpers2_20260210T174814.mat');
data10 = data10_t.rawData;
plot_digitizer = load("plot-data-cyl.csv");
plot_digitizer_X = plot_digitizer(:,1)';
plot_digitizer_Y = plot_digitizer(:,2)';

% Given data and constant
Patm = 102641.44; %Pa
Tatm = 294.65; %K
R = 287.05; % J/Kg K
row_atm = Patm/(R*Tatm); %Kg/m^3
q_10 = 60.7; %Pa
q_15 = 136.575; %Pa
mule = 1.716*10^-5*(Tatm/273.15)^(3/2)*((273.15+110.4)/(Tatm+110.4)); % Sutherland's
cyl_d = 0.0356; %m
conv = 248.84; %inH20 to Pa
dPatm = 0.025 * 3386.39; % Pa
dTatm = 0.5; % K
thetas = [-10,-5,0,5,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190];
theta_uncertantiy = 0.5;

% Calculate density uncertanity
Drow_atm_P = (Patm+dPatm)/(R*Tatm)-row_atm;
Drow_atm_T = Patm/(R*(Tatm+dTatm))-row_atm;
Drow_atm = abs(Drow_atm_T)+abs(Drow_atm_P);

% Calaculate free-stream velocity
mean_deg0_15 = mean(data15(:,4));
mean_deg0_10 = max(data10(:,4));
dv_15 = 15-sqrt((2*(mean_deg0_15+q_15))/row_atm);
dv_10 = 10-sqrt((2*(mean_deg0_10+q_10))/row_atm);

% Calculate Diferential preassures for each each velocity vs theta
DP_15 = (mean(data15(:,2:end)-mean(data15(:,1))))*conv;
DP_10 = (mean(data10(:,2:end)-mean(data10(:,1))))*conv;
dDP_15 = ((1.96*std(data15(:,2:end)-mean(data15(:,1))))/sqrt(25))*conv;
dDP_10 = ((1.96*std(data10(:,2:end)-mean(data10(:,1))))/sqrt(25))*conv;
figure(1)
errorbar(thetas,DP_15,dDP_15,'-o','LineWidth',.1,'MarkerSize',8,'MarkerFaceColor','r');
xlabel('Cylinder AoA (deg)');
ylabel('Differential Preassure (Pa)');
title('Measured Differential Preassure of Cylinder Surface Versus Angle of Attack at Wind Tunnel Air Speed of 15m/s','FontSize',14,'FontWeight','bold');
grid on;
figure(2)
errorbar(thetas,DP_10,dDP_10,'-o','LineWidth',.1,'MarkerSize',8,'MarkerFaceColor','r');
xlabel('Cylinder AoA (deg)');
ylabel('Differential Preassure (Pa)');
title('Measured Differential Preassure of Cylinder Surface Versus Angle of Attack at Wind Tunnel Air Speed of 10m/s','FontSize',14,'FontWeight','bold');
grid on;

% Calculate the coefficient of pressure for each velocity
theory_CP = (1-4*sin(deg2rad(thetas)).^2);
CP_15 = DP_15/max(DP_15);
CP_10 = DP_10/max(DP_10);

dCP_15_DP = ((DP_15+dDP_15)/max(DP_15))-CP_15;
dCP_15_dV = DP_15/(0.5*row_atm*((15+dv_15)^2))-CP_15;
dCP_15_drow = DP_15/(0.5*(row_atm+Drow_atm)*(15^2))-CP_15;
dCP_15 = abs(dCP_15_drow)+abs(dCP_15_dV)+abs(dCP_15_DP);

dCP_10_DP = ((DP_10+dDP_10)/max(DP_10))-CP_10;
dCP_10_dV = DP_10/(0.5*row_atm*((10+dv_10)^2))-CP_10;
dCP_10_drow = DP_10/(0.5*(row_atm+Drow_atm)*(10^2))-CP_10;
dCP_10 = abs(dCP_10_drow)+abs(dCP_10_dV)+abs(dCP_10_DP);

figure(3)
errorbar(thetas(:,3:23),CP_15(:,3:23),dCP_15(:,3:23),'-o','LineWidth',.1,'MarkerSize',8,'MarkerFaceColor','r'); hold on;
plot(thetas(:,3:23),theory_CP(:,3:23),'g','LineWidth',.05,'MarkerSize',5,'MarkerFaceColor','g'); hold off;
xlabel('Cylinder AoA (deg)');
ylabel('Coefficient of Preassure');
legend 'Measure CP' 'Theory CP';
title('Surface Coefficient of Preassure of Cylinder Versus Angle of Attack at Wind Tunnel Air Speed of 15m/s','FontSize',14,'FontWeight','bold');
grid on;

figure(4)
errorbar(thetas(:,3:23),CP_10(:,3:23),dCP_10(:,3:23),'-o','LineWidth',.1,'MarkerSize',8,'MarkerFaceColor','r'); hold on;
plot(thetas(:,3:23),theory_CP(:,3:23),'g','LineWidth',.05,'MarkerSize',5,'MarkerFaceColor','g'); hold off;
xlabel('Cylinder AoA (deg)');
ylabel('Coefficient of Preassure');
legend 'Measure CP' 'Theory CP';
title('Surface Coefficient of Preassure of Cylinder Versus Angle of Attack at Wind Tunnel Air Speed of 10m/s','FontSize',14,'FontWeight','bold');
grid on;

% Calculate Reynolds Number
Re_15 = (row_atm*(sqrt((2*max(data15))/row_atm))*cyl_d)/mule;
Re_10 = (row_atm*(sqrt((2*max(data10))/row_atm))*cyl_d)/mule;

% Calculate Cd
L_CP_15 = length(CP_15);
INT_theta = linspace(0,pi,L_CP_15);
integrand_15 = CP_15.*cos(INT_theta);
Cd_15 = trapz(INT_theta,integrand_15);
integrand_10 = CP_10.*cos(INT_theta);
Cd_10 = trapz(INT_theta,integrand_10);

dINT_theta = linspace(0.5,pi,L_CP_15);
dtheta_integrand_15 = CP_15.*cos(dINT_theta);
dCd_15_theta = trapz(dINT_theta,dtheta_integrand_15)-Cd_15;
dtheta_integrand_10 = CP_10.*cos(dINT_theta);
dCd_10_theta = trapz(dINT_theta,dtheta_integrand_10)-Cd_10;

dcp_integrand_15 = (CP_15+dCP_15).*cos(INT_theta);
dCd_15_cp = trapz(INT_theta,dcp_integrand_15)-Cd_15;
dcp_integrand_10 = (CP_10+dCP_10).*cos(INT_theta);
dCd_10_cp = trapz(INT_theta,dcp_integrand_10)-Cd_10;

dCd_15 = abs(dCd_15_theta)+abs(dCd_15_cp);
dCd_10 = abs(dCd_10_theta)+abs(dCd_10_cp);

% Calculate 2 dimensional Drag
Drag_15 = 0.5*row_atm*(15^2)*Cd_15*cyl_d;
Drag_10 = 0.5*row_atm*(10^2)*Cd_10*cyl_d;

figure(5)
errorbar(max(Re_15),Cd_15,dCd_15,'-o','LineWidth',.3,'MarkerSize',8,'MarkerFaceColor','r'); hold on;
errorbar(max(Re_10),Cd_10,dCd_10,'-o','LineWidth',.3,'MarkerSize',8,'MarkerFaceColor','c');
plot(plot_digitizer_X,plot_digitizer_Y,'g','LineWidth',.05,'MarkerSize',5,'MarkerFaceColor','g'); hold off;
xlabel('Reynolds Number');
ylabel('Coefficient of drag');
set(gca, 'XScale', 'log'); 
xlim([1e-1 1e7]);
xticks([1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]);
set(gca, 'YScale', 'log'); 
ylim([0.1 100]);
yticks([0.1 1 10 100]);
axis tight;
legend '15m/s Coefficient of Drag' '10m/s Coefficient of Drag' 'Cylinder Theorical Coefficient of Drag Curve';
title('Total Coefficeint of Drag of cylinder at Wind Tunnel Air Speeds of 10m/s and 15m/s Compared to Theorical Data','FontSize',14,'FontWeight','bold');
grid on;
