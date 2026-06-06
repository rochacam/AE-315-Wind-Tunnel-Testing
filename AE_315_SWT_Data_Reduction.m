% Code Name: SWT Aerodynamics Lab Data Reductor
% Code Description: Reduces the data collected in the SWT lab and plots
% necessary graphs and tables
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE315 - Section 29DB
% Date: 03/31/2026 until 04/05/2026
% Worked With: N/a

% Iniciation
clear; clc; close all;
raw_dta = readmatrix("DB29TEst2.csv");
Preassure_dta = raw_dta(:,2:17);

% Given data and constant
conv_length = 0.0254; %in to m
conv_preassure = 6894.76; %PSI to Pa
Patm = 101760; %Pa
Tatm = 294.65; %K
T0 = 296.15; %K
gamma = 1.4;
R = 287; % J/Kg K
row_atm = Patm/(R*Tatm); %Kg/m^3
A1_over_At = 3.27;
depth = 1*conv_length; %tunnel depth
A16_over_At = 1.7533;
At = 0.75*(conv_length^2);
Vatm = 15; %m/s
dPatm = 0.025 * 3386.39; % Pa
dTatm = 0.5; % K
port_locations = [1.5, 3, 3.75, 4.5, 5.25, 6, 6.75, 7.5, 8.25, 9, 9.75, 10.5, 11.25, 12, 12.75, 13.5].*conv_length;

% Finding Range of Steady State
figure(1);
plot(raw_dta(:,1),Preassure_dta(:,1),'-r','LineWidth',2); axis equal; grid on;
xlabel('Frame'); ylabel('Gauge Preassure(Pa)');
title('First Preassure Scanner Collected Data per Frame','FontSize',14,'FontWeight','bold');

% Ploting Average Preassure of each port vs x-location of each port
Steady_Preassure_dta = Preassure_dta(2750:3750,:).*conv_preassure;
Pabs = Steady_Preassure_dta+Patm;
Ps = mean(Pabs);
dSteadyP = ((1.96*std(Steady_Preassure_dta)-mean(Steady_Preassure_dta))/sqrt(2001));
DPs = abs(zeros(1,length(Ps))+dPatm)+abs(dSteadyP);
figure(2);
errorbar(port_locations,Ps,DPs,'vertical','-r','LineWidth',2); grid on;
xlabel('Port Location (m)'); ylabel('Steady Static Preassure (Pa)');
title('Static Preassure Distribution per Preassure Scanner Port','FontSize',14,'FontWeight','bold');

% Compressible Flow Calculator Properties for station 1 and 16
P1_over_Po_1 = 0.97753557; % Data from calculator
T1_over_T0 = 0.99352942; % Data from calculator
M1_sub =  0.18045380; % Data from calculator
Po = Ps(:,1)/P1_over_Po_1
NPR = Po/Patm
T1 = T1_over_T0*T0;

P16_over_P0 =  0.11905949; % Data from calculator
T16_over_T0 =  0.54441813; % Data from calculator
M16_sup =  2.04551151; % Data from calculator
P0_16 = Ps(:,16)/P16_over_P0;
Pe_over_P0 = Patm/P0_16;

% Shock Strength (Shock between port 11 and 12)
P_before_NS = Ps(:,11);
P_after_NS = Ps(:,12);
SS = P_after_NS/P_before_NS

% Find Mass flow rate
m_dot = (At*(Po/sqrt(T0)))*(sqrt(gamma/R)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))))

% Find Ideal Thrust of Nozzle
V_inlet = M1_sub*sqrt(gamma*R*T1);
v_exit = sqrt(((2*gamma*R*T0)/(gamma-1))*(1-(Pe_over_P0)^((gamma-1)/gamma))+V_inlet^2);

T = m_dot*(v_exit-V_inlet)
