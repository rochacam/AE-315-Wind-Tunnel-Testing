% Code Name: INT Aerodynamics Lab Data Reductor
% Code Description: Reduces the data collected in the INT lab and plots
% necessary graphs and tables
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE315 - Section 29DB
% Date: 02/02/2026
% Worked With: N/a

% Iniciation
clear; clc; close all;
data1 = load("0Hz_20260127T180458.mat");
DP1 = data1.rawData;


% Given data and constant
Patm = 102946.215; %Pa
Tatm = 295.15; %K
R = 287.05; % J/Kg K
row_atm = Patm/(R*Tatm); %Kg/m^3
q = 137.81; %Pa
mule = 1.716*10^-5*(Tatm/273.15)^(3/2)*((273.15+110.4)/(Tatm+110.4)); % Sutherland's
t = 2.093;
conv = 248.84; %inH20 to Pa

% Uncertanties
actual_dP = mean(DP1-DP1(:,1))*conv;
dPatm = 0.025 * 3386.39;
dTatm = 0.5;

% Data Separation
P0 = DP1(:,1);
P4 = DP1(:,2);
P8 = DP1(:,3);
P12 = DP1(:,4);
P18 = DP1(:,5);
group = {P4,P8,P12,P18};
names = [4,8,12,18];
mean_P0 = mean(P0);
std_P0 = std(P0);

% Uncertanty Propagation
mean_q = zeros(1,4);
std_q = zeros(1,4);
sem_q = zeros(1,4);
deltaP = zeros(1,4);
air_speed = zeros(1,4);
dV = zeros(1,4);
Re = zeros(1,4);

for i = 1:4
    m = mean(group{i})-mean_P0;
    s = std(group{i});
    N = length(group{i});
    SEM = s/sqrt(N);
    mean_q(i) = m*conv;
    std_q(i) = s*conv;
    SEM_q(i) = SEM*conv;
    deltaP(i) = t*SEM*conv;
    air_speed(i) = sqrt(2*mean_q(i)*R*Tatm/Patm);
    V = air_speed(i);
    P_q = mean_q(i);
    dP_q = deltaP(i);

    dV_P_q = sqrt(2*(P_q+dP_q)*R*Tatm/Patm)-V;
    dV_P = sqrt(2*R*Tatm*P_q/(Patm+dPatm))-V;
    dV_T = sqrt(2*R*P_q*(Tatm-dTatm)/Patm)-V;
    dV(i) = abs(dV_T)+abs(dV_P_q)+abs(dV_P);

    Re(i) = (row_atm*air_speed(i)/mule);
end

Table_1 = table(names(:), mean_q(:), deltaP(:), air_speed(:), dV(:), Re(:),'VariableNames',{'Frqeuency in Hz', 'Mean Dynamic Preassure', 'Delta Pressure','airspeed','Delta velocity','Reynolds Number'});
disp(Table_1);

figure;
errorbar(names,air_speed,dV,'-o','LineWidth',.1,'MarkerSize',8,'MarkerFaceColor','r');
xlabel('Fan Frequency (Hz)');
ylabel('Air Speed (m/s)');
title('Wind Tunnel Air Speed vs Fan Frequency','FontSize',14,'FontWeight','bold');
grid on;

% Calculations for the 40 points 0Hz and 20Hz
% Traverse
load("0HZ_20260127T174219.mat")
tare = rawData;
load("20HZ_20260127T175122.mat")
data20 = rawData;

% Calculate average q_infinity
q_infinity = mean((data20'-mean(tare')))*conv;

% Uncertanty of q_infinity with 95% confidence
delta_q = (1.96*(std(data20'-mean(tare')))/sqrt(2000))*conv;
heights = linspace(0,8,40)'*0.0254;

% Velocity at each q_infinity calculation
velocity = sqrt(2*q_infinity/row_atm);
velocity_uncertanity = zeros(size(velocity));
for j =1:length(velocity)
    VL = velocity(j);
    q_20 = q_infinity(j);
    DP_20 = delta_q(j);
    DVL_Pq = sqrt(2*(q_20+DP_20)/row_atm)-VL;
    DVL_Patm = sqrt(2*q_20*R*Tatm/(Patm+dPatm))-VL;
    DVL_Tatm = sqrt(2*q_20*R*(Tatm+dTatm)/Patm)-VL;
    velocity_uncertanity(j) = abs(DVL_Tatm)+abs(DVL_Pq)+abs(DVL_Patm);
end

I = zeros(size(velocity));

for k = 1:length(velocity)
    VLO = velocity(k);
    I(k) = (std(velocity)/VLO)*100;
end
mean_I = mean(I);
std_I = std(I);

% Plot Velocity Vs height
figure(3);
errorbar(velocity,heights,velocity_uncertanity,'horizontal','o-','LineWidth',1.5)
xlabel('Mean Velocity (m/s)');
ylabel('Heights (m)');
title('Velocity profile in the wind tunnel','FontSize',14,'FontWeight','bold');
ylim([-0.2 0.25]);

% Plot Turbulance intenisty vs height
figure(4);
plot(heights,I,'o-')
yline(mean_I,'--r','Mean Turbulance Intencity')
xlabel('Heights (m)');
ylabel('Turbulance Intensity (%)');
title('Turbulance Intensity versus tunnel height','FontSize',14,'FontWeight','bold');
grid on;


