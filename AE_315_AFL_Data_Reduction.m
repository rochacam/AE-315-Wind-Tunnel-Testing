% Code Name: AFL Aerodynamics Lab Data Reductor
% Code Description: Reduces the data collected in the AFL lab and plots
% necessary graphs and tables
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE315 - Section 29DB
% Date: 03/03/2026 until
% Worked With: N/a

% Iniciation
clear; clc; close all;
conv = 248.84; %inH20 to Pa
trav_tare_dta = load('15deg_traverse_tare_20260224T173614.csv');
trav_neg10_dta = (load('neg10deg_traverse_15ms_20260224T175407.csv')-trav_tare_dta).*conv;
trav_neg5_dta = (load('neg5deg_traverse_15ms_20260224T180400.csv')-trav_tare_dta).*conv;
trav_0_dta = (load('zero_deg_traverse_15ms_20260224T181255.csv')-trav_tare_dta).*conv;
trav_5_dta = (load('five_deg_traverse_15ms_20260224T182138.csv')-trav_tare_dta).*conv;
trav_10_dta = (load('ten_deg_traverse_15ms_20260224T183038.csv')-trav_tare_dta).*conv;
trav_15_dta = (load('fifteen_deg_traverse_15ms_20260224T183846.csv')-trav_tare_dta).*conv;

% Given data and constant
Patm = 102608; %Pa
Tatm = 294.65; %K
R = 287; % J/Kg K
row_atm = Patm/(R*Tatm); %Kg/m^3
q = 136.46; %Pa
mule = 1.716*10^-5*(Tatm/273.15)^(3/2)*((273.15+110.4)/(Tatm+110.4)); % Sutherland's
Vatm = 15; %m/s
conv = 248.84; %inH20 to Pa
dPatm = 0.025 * 3386.39; % Pa
dTatm = 0.5; % K
thetas = [-10,-5,0,5,10,15];
theta_uncertantiy = 0.5;

% Setup Empty Running Totals
Cl_i = zeros(1,1);
Cmle_i = zeros(1,1);
Cmcby4_i = zeros(1,1);
CL_exp_i = zeros(1,1);
dCL_exp_i = zeros(1,1);
Drag_i = zeros(1,1);
d_Drag_i = zeros(1,1);
CD_i = zeros(1,1);
d_CD_i = zeros(1,1);
CmLE_exp_i = zeros(1,1);
Cmc4_exp_i = zeros(1,1);
CD_Panel_i = zeros(1,1);

% Airfoil preassure tube locations
x_loc_top = [0.0000,0.1500,0.3000,0.6000,1.2000,1.7998,2.4001,3.0000,4.2000,5.4000]./6;
x_loc_bottom = [0.0000,0.2400,0.3600,0.6000,1.2000,1.8000,3.000,4.2000,5.1500]./6;
x_over_c_top = [0.0000,0.1500,0.3000,0.6000,1.2000,1.7998,2.4001,3.0000,4.2000,5.4000,6]./6;
x_over_c_bottom =  [0.0000,0.2400,0.3600,0.6000,1.2000,1.8000,3.000,4.2000,5.1500,6]./6;
x_over_c = sortrows([x_over_c_top(2:10),x_over_c_bottom]',"ascend");
y_trav = linspace(-10.16,10.16,40);

% Wake Profile at each AoA
trav_neg10 = mean(trav_neg10_dta,2);
trav_neg5 = mean(trav_neg5_dta,2);
trav_0 = mean(trav_0_dta,2);
trav_5 = mean(trav_5_dta,2);
trav_10 = mean(trav_10_dta,2);
trav_15 = mean(trav_15_dta,2);

v_wake_neg10 = sqrt(2*(trav_neg10)/row_atm);
v_wake_neg5 = sqrt(2*(trav_neg5)/row_atm);
v_wake_0 = sqrt(2*(trav_0)/row_atm);
v_wake_5 = sqrt(2*(trav_5)/row_atm);
v_wake_10 = sqrt(2*(trav_10)/row_atm);
v_wake_15 = sqrt(2*(trav_15)/row_atm);
v_wake_all = [v_wake_neg10,v_wake_neg5,v_wake_0,v_wake_5,v_wake_10,v_wake_15];

% Velocity wake Uncertanty
Drow_atm_P = (Patm+dPatm)/(R*Tatm)-row_atm;
Drow_atm_T = Patm/(R*(Tatm+dTatm))-row_atm;
Drow_atm = abs(Drow_atm_T)+abs(Drow_atm_P);

dtrav_neg10 = ((1.96*std(trav_neg10)-mean(trav_neg10))/sqrt(5000));
dtrav_neg5 = ((1.96*std(trav_neg5)-mean(trav_neg5))/sqrt(5000));
dtrav_0 = ((1.96*std(trav_0)-mean(trav_0))/sqrt(5000));
dtrav_5 = ((1.96*std(trav_5)-mean(trav_5))/sqrt(5000));
dtrav_10 = ((1.96*std(trav_10)-mean(trav_10))/sqrt(5000));
dtrav_15 = ((1.96*std(trav_15)-mean(trav_15))/sqrt(5000));

dv_wake_neg10 = abs((sqrt(2*(trav_neg10+dtrav_neg10)/row_atm))-v_wake_neg10)+abs(sqrt(2*(trav_neg10)/(row_atm+Drow_atm))-v_wake_neg10);
dv_wake_neg5 = abs((sqrt(2*(trav_neg5+dtrav_neg5)/row_atm))-v_wake_neg5)+abs(sqrt(2*(trav_neg5)/(row_atm+Drow_atm))-v_wake_neg5);
dv_wake_0 = abs((sqrt(2*(trav_0+dtrav_0)/row_atm))-v_wake_0)+abs(sqrt(2*(trav_0)/(row_atm+Drow_atm))-v_wake_0);
dv_wake_5 = abs((sqrt(2*(trav_5+dtrav_5)/row_atm))-v_wake_5)+abs(sqrt(2*(trav_5)/(row_atm+Drow_atm))-v_wake_5);
dv_wake_10 = abs((sqrt(2*(trav_10+dtrav_10)/row_atm))-v_wake_10)+abs(sqrt(2*(trav_10)/(row_atm+Drow_atm))-v_wake_10);
dv_wake_15 = abs((sqrt(2*(trav_15+dtrav_15)/row_atm))-v_wake_15)+abs(sqrt(2*(trav_15)/(row_atm+Drow_atm))-v_wake_15);
dv_wake_all = [dv_wake_neg10, dv_wake_neg5, dv_wake_0, dv_wake_5, dv_wake_10, dv_wake_15];

% Velocity Wake Plots
errorbar(v_wake_15,y_trav,dv_wake_15,'horizontal','g'); hold on;
errorbar(v_wake_10,y_trav,dv_wake_10,'horizontal','b');
errorbar(v_wake_5,y_trav,dv_wake_5,'horizontal','c');
errorbar(v_wake_0,y_trav,dv_wake_0,'horizontal','y');
errorbar(v_wake_neg5,y_trav,dv_wake_neg5,'horizontal','m');
errorbar(v_wake_neg10,y_trav,dv_wake_neg10,'horizontal','r'); hold off;
xlabel('Wake Velocity (m/s)');
ylabel('Wind Tunnel Height (cm)');
legend 'Wake Velocity at 15 deg AoA' 'Wake Velocity at 10 deg AoA' 'Wake Velocity at 5 deg AoA' 'Wake Velocity at 0 deg AoA' 'Wake Velocity at -5 deg AoA' 'Wake Velocity at -10 deg AoA';
title('Wake Velocity behind NACA 4412 airfoil versus Angle of Attack at Wind Tunnel Air Speed of 15m/s','FontSize',14,'FontWeight','bold');
grid on;

% Wake points Begin/end for each AOA
y_neg10 = linspace((-5.47077),(4.42872),40);
y_neg5 = linspace((0.260513),(2.34462),40);
y_0 = linspace((-1.30256),(0.781538),40);
y_5 = linspace((-2.86564),(0.260513),40);
y_10 = linspace((-3.38667),(1.302566),40);
y_15 = linspace((-4.42872),(4.42872),40);
y_start_end = [y_neg10', y_neg5', y_0', y_5', y_10', y_15'];

% Get Free-Stream velocity and dynamic preassure at each AoA
trav_neg10(10:28) = [];
trav_neg5(20:25) = [];
trav_0(18:22) = [];
trav_5(15:21) = [];
trav_10(14:23) = [];
trav_15(12:29) = [];
v_inf_neg10 = mean(sqrt(2*(trav_neg10)/row_atm));
v_inf_neg5 = mean(sqrt(2*(trav_neg5)/row_atm));
v_inf_0 = mean(sqrt(2*(trav_0)/row_atm));
v_inf_5 = mean(sqrt(2*(trav_5)/row_atm));
v_inf_10 = mean(sqrt(2*(trav_10)/row_atm));
v_inf_15 = mean(sqrt(2*(trav_15)/row_atm));
v_inf_all = [v_inf_neg10, v_inf_neg5, v_inf_0, v_inf_5, v_inf_10, v_inf_15];
q_neg10 = 0.5 * row_atm * v_inf_neg10^2;
q_neg5 = 0.5 * row_atm * v_inf_neg5^2;
q_0 = 0.5 * row_atm * v_inf_0^2;
q_5 = 0.5 * row_atm * v_inf_5^2;
q_10 = 0.5 * row_atm * v_inf_10^2;
q_15 = 0.5 * row_atm * v_inf_15^2;
q_all = [q_neg10,q_neg5,q_0,q_5,q_10,q_15];

% Calculate uncertantiy for dynamic preassures
dtrav_neg10 = ((1.96*std(trav_neg10)-mean(trav_neg10))/sqrt(40));
dtrav_neg5 = ((1.96*std(trav_neg5)-mean(trav_neg5))/sqrt(5000));
dtrav_0 = ((1.96*std(trav_0)-mean(trav_0))/sqrt(5000));
dtrav_5 = ((1.96*std(trav_5)-mean(trav_5))/sqrt(5000));
dtrav_10 = ((1.96*std(trav_10)-mean(trav_10))/sqrt(5000));
dtrav_15 = ((1.96*std(trav_15)-mean(trav_15))/sqrt(5000));

dv_inf_neg10 = abs((sqrt(2*(trav_neg10+dtrav_neg10)/row_atm))-v_inf_neg10)+abs(sqrt(2*(trav_neg10)/(row_atm+Drow_atm))-v_inf_neg10)-1;
dv_inf_neg5 = abs((sqrt(2*(trav_neg5+dtrav_neg5)/row_atm))-v_inf_neg5)+abs(sqrt(2*(trav_neg5)/(row_atm+Drow_atm))-v_inf_neg5);
dv_inf_0 = abs((sqrt(2*(trav_0+dtrav_0)/row_atm))-v_inf_0)+abs(sqrt(2*(trav_0)/(row_atm+Drow_atm))-v_inf_0);
dv_inf_5 = abs((sqrt(2*(trav_5+dtrav_5)/row_atm))-v_inf_5)+abs(sqrt(2*(trav_5)/(row_atm+Drow_atm))-v_inf_5);
dv_inf_10 = abs((sqrt(2*(trav_10+dtrav_10)/row_atm))-v_inf_10)+abs(sqrt(2*(trav_10)/(row_atm+Drow_atm))-v_inf_10);
dv_inf_15 = abs((sqrt(2*(trav_15+dtrav_15)/row_atm))-v_inf_15)+abs(sqrt(2*(trav_15)/(row_atm+Drow_atm))-v_inf_15);
dv_inf_all = [mean(dv_inf_neg10)', mean(dv_inf_neg5)', mean(dv_inf_0)', mean(dv_inf_5)', mean(dv_inf_10)', mean(dv_inf_15)'];

dq_neg10 = mean(dv_inf_neg10 + Drow_atm);
dq_neg5 = mean(dv_inf_neg5 + Drow_atm);
dq_0 = mean(dv_inf_0 + Drow_atm);
dq_5 = mean(dv_inf_5 + Drow_atm);
dq_10 = mean(dv_inf_10 + Drow_atm);
dq_15 = mean(dv_inf_15 + Drow_atm);
dq_all = [dq_neg10,dq_neg5,dq_0,dq_5,dq_10,dq_15];

% Calcuate 2-Dimensional Drag for each AoA
for k = 1:length(v_inf_all)
    Integrand = v_wake_all(:,k).*(v_inf_all(:,k)-v_wake_all(:,k));
    Drag = trapz(y_neg10,Integrand);
    D_Integrand = (v_wake_all(:,k)+dv_wake_all(:,k)).*((v_inf_all(:,k)+dv_inf_all(:,k))-(v_wake_all(:,k)+dv_wake_all(:,k)));
    d_Drag = trapz(y_start_end(:,k),D_Integrand);
    Cd = Drag/(q_all(:,k)*10.16);
    d_Cd = abs(((Drag+((d_Drag-Drag).*10^-1))/(q_all(:,k)*10.16))-Cd)+abs((Drag/((q_all(:,k)+dq_all(:,k))*10.16))-Cd);
    Drag_i = [Drag_i, Drag];
    d_Drag_i = [d_Drag_i, (d_Drag-Drag).*10^-1];
    CD_i = [CD_i, Cd];
    d_CD_i = [d_CD_i, d_Cd];
end

% Load and slice Preassure Scanner Data
ps_neg10_dta = readcell('neg10deg_pressure_scanner_15ms_AFL.csv');
ps_neg5_dta = readcell('neg5deg_pressure_scanner_15ms_AFL.csv');
ps_0_dta = readcell('zero_deg_pressure_scanner_15ms_AFL.csv');
ps_5_dta = readcell('five_deg_pressure_scanner_15ms_AFL.csv');
ps_10_dta = readcell('ten_deg_pressure_scanner_15ms_AFL.csv');
ps_15_dta = readcell('fifteen_deg_pressure_scanner_15ms_AFL.csv');

fs_ps_neg10 = cell2mat(ps_neg10_dta(2:end,2));
fs_ps_neg5 = cell2mat(ps_neg5_dta(2:end,2));
fs_ps_0 = cell2mat(ps_0_dta(2:end,2));
fs_ps_5 = cell2mat(ps_5_dta(2:end,2));
fs_ps_10 = cell2mat(ps_10_dta(2:end,2));
fs_ps_15 = cell2mat(ps_15_dta(2:end,2));

ports_ps_neg10 = mean(cell2mat(ps_neg10_dta(2:end,3:20)),1);
ports_ps_neg5 = mean(cell2mat(ps_neg5_dta(2:end,3:20)),1);
ports_ps_0 = mean(cell2mat(ps_0_dta(2:end,3:20)),1);
ports_ps_5 = mean(cell2mat(ps_5_dta(2:end,3:20)),1);
ports_ps_10 = mean(cell2mat(ps_10_dta(2:end,3:20)),1);
ports_ps_15 = mean(cell2mat(ps_15_dta(2:end,3:20)),1);

cp_neg10 = [(ports_ps_neg10.*conv)/q_neg10, 0];
cp_neg5 = [(ports_ps_neg5.*conv)/q_neg5, 0];
cp_0 = [(ports_ps_0.*conv)/q_0, 0];
cp_5 = [(ports_ps_5.*conv)/q_5, 0];
cp_10 = [(ports_ps_10.*conv)/q_10, 0];
cp_15 = [(ports_ps_15.*conv)/q_15, 0];
cp_aoa = [cp_neg10', cp_neg5', cp_0', cp_5', cp_10', cp_15'];

dcp_neg10 = [((ports_ps_neg10.*conv)/(q_neg10+dq_neg10))-cp_neg10(1:18), 0];
dcp_neg5 = [((ports_ps_neg5.*conv)/(q_neg5+dq_neg5))-cp_neg5(1:18), 0];
dcp_0 = [((ports_ps_0.*conv)/(q_0+dq_0))-cp_0(1:18), 0];
dcp_5 = [((ports_ps_5.*conv)/(q_5+dq_5))-cp_5(1:18), 0];
dcp_10 = [((ports_ps_10.*conv)/(q_10+dq_10))-cp_10(1:18), 0];
dcp_15 = [((ports_ps_15.*conv)/(q_15+dq_15))-cp_15(1:18), 0];
dcp_aoa = [dcp_neg10', dcp_neg5', dcp_0', dcp_5', dcp_10', dcp_15'];

% Panel Method code
for i = 1:length(thetas)
    af_no=4412;
    no_pls=201;
    alpha=thetas(1,i);
    V_inf=1;
    [x,z,xu,xl]=NACA4(af_no,no_pls);
    no_pls=length(x);
    [sth,cth]=setUpPanel1(x,z);
    [coeffs,beta,r]=calCoeffs1(x,z,sth,cth,alpha,V_inf);
    [Vt]=calVt(sth,cth,r,beta,V_inf,alpha,coeffs);
    Cp=1-(Vt./V_inf).^2;
    Cpl=Cp(1:length(xl));
    Cpu=Cp(length(xl)+1:end);
    xu=[xu xl(1)];
    Cpu=[Cpu Cpl(1)];
    xl=[xl xu(1)];
    Cpl=[Cpl Cpu(1)];
    Cl=simps(fliplr(xl),fliplr(Cpl))-simps(xu,Cpu);
    CmLE=simps(xu,Cpu.*xu)-simps(fliplr(xl),fliplr(Cpl.*xl));
    xCp=-CmLE/Cl;
    CmcBy4=Cl/4+CmLE;

    Cl_i = [Cl_i, Cl];
    Cmle_i = [Cmle_i, CmLE];
    Cmcby4_i = [Cmcby4_i, CmcBy4];

    cp_i_upper = cp_aoa([1:10 19],i);
    cp_i_lower = cp_aoa([1 11:18 19],i);
    dcp_i_upper = dcp_aoa([1:10 19],i);
    dcp_i_lower = dcp_aoa([1 11:18 19],i);

    top_iternpolation = interp1(x_over_c_top,cp_i_upper,0.96,'linear','extrap');
    bottom_iternpolation = interp1(x_over_c_bottom,cp_i_lower,0.96,"linear",'extrap');
    TE_cp_i = (top_iternpolation+bottom_iternpolation)/2;
    d_top_iternpolation = interp1(x_over_c_top,dcp_i_upper,0.96,'linear','extrap');
    d_bottom_iternpolation = interp1(x_over_c_bottom,dcp_i_lower,0.96,"linear",'extrap');
    dTE_cp_i = (d_top_iternpolation+d_bottom_iternpolation)/2;

    cp_lower = [cp_i_lower(1:9,:); TE_cp_i];
    cp_upper = [cp_i_upper(1:10,:); TE_cp_i];
    dcp_lower = [dcp_i_lower(1:9,:); dTE_cp_i];
    dcp_upper = [dcp_i_upper(1:10,:); dTE_cp_i];

    CL_exp = (1/6)*(trapz(x_over_c_bottom.*6,cp_lower)-trapz(x_over_c_top.*6,cp_upper));
    CL_exp_i = [CL_exp_i, CL_exp];
    dCL_exp = ((1/6)*(trapz(x_over_c_bottom.*6,(cp_lower+dcp_lower)*1.02)-trapz(x_over_c_top.*6,(cp_upper+dcp_upper)*1.02)))-CL_exp;
    dCL_exp_i = [dCL_exp_i, dCL_exp];

    CD_Panel = (CL_exp^2)/(pi*12*0.76);
    CD_Panel_i = [CD_Panel_i, CD_Panel];

    CmLE_exp = -(1/6^2)*(trapz(x_over_c_bottom.*6,(cp_lower*(x_loc_bottom.*6)))-trapz(x_over_c_top.*6,(cp_upper*(x_loc_top(2:end).*6))));
    CmLE_exp_i = [CmLE_exp_i, mean(CmLE_exp(2:end))];

    Cmc4_exp = -(1/6^2)*(trapz(x_over_c_bottom.*6,(cp_lower*((x_loc_bottom*6)-1.5)))-trapz(x_over_c_top.*6,(cp_upper*((x_loc_top(2:end)*6)-1.5))));
    Cmc4_exp_i = [Cmc4_exp_i, mean(Cmc4_exp(2:end))];

    figure(1+i);
    plot(xu,Cpu,'-r',xl,Cpl,'-b','LineWidth',2); hold on;
    errorbar(x_over_c_bottom,cp_lower,dcp_lower.*2,'vertical','m','LineWidth',2,"MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]);
    errorbar(x_over_c_top,cp_upper,dcp_upper.*2,'vertical','y','LineWidth',2,"MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]); hold off;
    hold on
    dynamic_title = sprintf('Surface Coefficient of Pressure distribution for NACA 4412 airfoil at an AoA of %d', thetas(1,i));
    title(dynamic_title,'FontSize',14,'FontWeight','bold');
    legend('Panel Method Upper','Panel Method Lower','EXP Lower','EXP Upper')
    xlabel('$x/c$','Interpreter','latex')
    ylabel('$C_p$','Interpreter','latex')
    set(gca,'YDir','reverse'); grid on;
end
% Linear Regression for Zero-Lift AoA
fit = polyfit(thetas,Cl_i(2:end),1);
m = fit(1);
b = fit(2);
AoA_zl = -b/m;

figure(8);
plot(x,z,'-r','LineWidth',2)
axis equal
grid;
xlabel('$x/c$','Interpreter','latex')
ylabel('$z/c$','Interpreter','latex')
figure(9)
plot(thetas,Cl_i(2:end),'-r','LineWidth',2); hold on;
errorbar(thetas,CL_exp_i(2:end),dCL_exp_i(2:end),'vertical','m','LineWidth',2);
scatter(AoA_zl,0,'g','filled'); hold off;
ylabel('$C_l$','Interpreter','latex');
xlabel('Angle of Attack (deg)');
legend ('Theorical Coefficient of Lift','Experimental Coefficient of Lift','Zero-Lift Angle of Attack');
title('Coefficent of Lift versus Angle of Attack of NACA 4412 Airfoil','FontSize',14,'FontWeight','bold');
grid on;
figure(10)
errorbar(thetas,CD_i(2:end),d_CD_i(2:end),'vertical','r','LineWidth',2); hold on;
yline(0,'g','LineWidth',2); hold off;
ylabel('$C_d$','Interpreter','latex');
xlabel('Angle of Attack (deg)');
legend ('Experimental Coefficient of Drag','Theorical Coefficient of Drag');
title('Coefficent of Drag versus Angle of Attack of NACA 4412 Airfoil','FontSize',14,'FontWeight','bold');
grid on;
figure(11)
plot(thetas,Cmle_i(2:end),'r','LineWidth',2,'Marker','o');hold on;
plot(thetas,CmLE_exp_i(2:end),'m','LineWidth',2,'Marker','o');
xlabel('Angle of Attack (deg)');
ylabel('Coefficicient of Moment about Leading Edge');
legend ('Theorical Cm','Experimental Cm');
title('Coefficent of Moment at Leading Edge versus Angle of Attack of NACA 4412 Airfoil','FontSize',14,'FontWeight','bold');
grid on; hold off;
figure(12)
plot(thetas,Cmcby4_i(2:end),'r','LineWidth',2,'Marker','o');hold on;
plot(thetas,Cmc4_exp_i(2:end),'m','LineWidth',2,'Marker','o');
xlabel('Angle of Attack (deg)');
ylabel('Coefficicient of Moment about Quarter Chord');
legend ('Theorical Cm','Experimental Cm');
title('Coefficent of Moment at Quarter Chord versus Angle of Attack of NACA 4412 Airfoil','FontSize',14,'FontWeight','bold');
grid on; hold off;
figure(13)
plot(thetas,CL_exp_i(2:end)./CD_i(2:end),'r','LineWidth',2,'Marker','o');
xlabel('Angle of Attack (deg)');
ylabel('$C_l$ / $C_d$','Interpreter','latex');
legend ('Experimental CL/CD ');
title('Coefficent of Lift over Coefficient of Drag versus Angle of Attack of NACA 4412 Airfoil','FontSize',14,'FontWeight','bold');
grid on;
figure(14)
plot(CD_i(2:end),CL_exp_i(2:end),'r','LineWidth',2,'Marker','o'); hold on;
plot(CD_Panel_i(2:end),CL_exp_i(2:end),'g','LineWidth',2,'Marker','o'); hold off
xlabel('$C_d$','Interpreter','latex');
ylabel('$C_l$','Interpreter','latex');
legend ('Experimental Drag Polar','Panel Method Drag Polar');
title('Drag Polar of NACA 4412 Airfoil','FontSize',14,'FontWeight','bold');
grid on;

%% Functions
function [x,z,x_u,x_l,z_u,z_l]=NACA4(af_no,no_pls)

% Preallocate arrays.
x_c=linspace(0,pi,floor(no_pls/2));
x_c=(1-cos(x_c))./2;
z_c=linspace(0,pi,floor(no_pls/2));

% Extract airfoil data from NACA airfoil designation.
t=(af_no-floor(af_no/100)*100)/100; % maximum thickness in tenth's of chord
m=floor(af_no/1000)/100; %maximum thickness in tenth's of chord
p=(floor(af_no/100)-floor(af_no/1000)*10)/10; %position of maximum thickness in tenth's of chord

% Create airfoil thickness distribution.
zthick=(t/0.2).*(.2969.*sqrt(x_c)-.126.*x_c-.35160.*x_c.^2+...
    .2843.*x_c.^3-.1036.*x_c.^4);

% Create airfoil camberline.
for i=1:1:length(x_c)
    if x_c(i)<=p
        z_c(i)=m*(1/p)^2.*(2.*p.*x_c(i)-x_c(i).^2);
    else
        z_c(i)=m*(1/(p-1))^2.*((1-2*p)+2.*p.*x_c(i)-x_c(i).^2);
    end
end

% Create airfoil upper (suction side) and lower (pressure side) surfaces by
% adding thickness distribution to camberline.
z_u=zthick+z_c;
z_l=z_c-zthick;

% Re-arrange coordinates to enable easy plotting.
z_l=fliplr(z_l(2:end));
x_l=fliplr(x_c(2:end));

z_u=z_u(1:end-1);
x_u=x_c(1:end-1);

if m==0
    z_u(1)=0;
end

z=[z_l z_u]';
x=[x_l x_u]';

end

function [sth,cth,len]=setUpPanel1(x,z)

N=length(x);

len=zeros(1,N);
sth=len;
cth=len;

for i=1:1:N
    if i==length(x)
        mm=1;
    else
        mm=i+1;
    end
    len(i)=sqrt((x(mm)-x(i))^2+(z(mm)-z(i))^2);
    sth(i)=(z(mm)-z(i))/len(i);
    cth(i)=(x(mm)-x(i))/len(i);
end

end

function [coeffs,beta,r]=calCoeffs1(x,z,...
    sth,cth,alpha,V_inf)

A=zeros(length(x)+1,length(x)+1);
b=zeros(1,length(x)+1);

[r,beta]=calRBeta(x,z);

alpha=alpha/180*pi;

for i=1:1:length(x)
    sum=0;
    for j=1:1:length(x)
        if j==length(x)
            mm=1;
        else
            mm=j+1;
        end

        sd=sinDiff(sth,cth,i,j);
        cd=cosDiff(sth,cth,i,j);

        A(i,j)=sd*log(r(i,mm)/r(i,j))+cd*beta(i,j);
        sum=sum+cd*log(r(i,mm)/r(i,j))-sd*beta(i,j);


        A(i,j)=A(i,j)/(2*pi);
    end
    A(i,length(x)+1)=sum/(2*pi);

    s1=(cos(alpha)*sth(i)-sin(alpha)*cth(i));
    b(i)=V_inf*s1;
end

tt=[1 length(x)];

for j=1:1:length(x)
    if j==length(x)
        mm=1;
    else
        mm=j+1;
    end
    sum=0;
    for k=1:1:length(tt)
        sd=sinDiff(sth,cth,tt(k),j);
        cd=cosDiff(sth,cth,tt(k),j);

        sum=sum+(sd*beta(tt(k),j)-...
            cd*log(r(tt(k),mm)/r(tt(k),j)));
    end
    A(length(x)+1,j)=sum/(2*pi);
end


sum2=0;
for k=1:1:length(tt)
    sum1=0;
    for j=1:1:length(x)
        if j==length(x)
            mm=1;
        else
            mm=j+1;
        end
        sd=sinDiff(sth,cth,tt(k),j);
        cd=cosDiff(sth,cth,tt(k),j);

        sum1=sum1+(sd*log(r(tt(k),mm)/r(tt(k),j))+...
            cd*beta(tt(k),j));
    end
    sum2=sum2+sum1;
end

A(length(x)+1,length(x)+1)=sum2/(2*pi);
c1=(sin(alpha)*sth(1)+cos(alpha)*cth(1));
c2=(sin(alpha)*sth(length(x))+cos(alpha)*cth(length(x)));
b(length(x)+1)=-V_inf*c1-V_inf*c2;

coeffs=A\b';

end

function [r,beta]=calRBeta(x,z)

r=zeros(length(x),length(x));
beta=r;
for i=1:1:length(x)
    if i==length(x)
        xcntr=(x(1)+x(i))/2;
        zcntr=(z(1)+z(i))/2;
    else
        xcntr=(x(i+1)+x(i))/2;
        zcntr=(z(i+1)+z(i))/2;
    end

    for j=1:1:length(x)
        r(i,j)=sqrt((xcntr-x(j))^2+(zcntr-z(j))^2);

        if i~=j
            if j==length(x)
                k=1;
            else
                k=j+1;
            end
            temp1=(zcntr-z(k))*(xcntr-x(j))-...
                (xcntr-x(k))*(zcntr-z(j));
            temp2=(xcntr-x(k))*(xcntr-x(j))+...
                (zcntr-z(k))*(zcntr-z(j));
            beta(i,j)=atan2(temp1,temp2);
        else
            beta(i,j)=pi;
        end

    end
end
end

function [sinD]=sinDiff(sin_theta,cos_theta,i,j)

sinD=(cos_theta(j)*sin_theta(i)-sin_theta(j)*cos_theta(i));

end

function [cosD]=cosDiff(sin_theta,cos_theta,i,j)

cosD=(sin_theta(j)*sin_theta(i)+cos_theta(j)*cos_theta(i));

end

function [Vt]=calVt(sth,cth,r,beta,V_inf,alpha,coeffs)


q=coeffs(1:end-1);
gamma=coeffs(end);

alpha=alpha/180*pi;

N=size(beta,1);
Vt=zeros(1,N);
for i=1:1:N
    sum1=0;
    sum2=0;
    for j=1:1:N
        if j==N
            mm=1;
        else
            mm=j+1;
        end

        sd=sinDiff(sth,cth,i,j);
        cd=cosDiff(sth,cth,i,j);

        sum1=sum1+q(j)/(2*pi)*(sd*beta(i,j)-cd*log(r(i,mm)/r(i,j)));
        sum2=sum2+gamma/(2*pi)*(sd*log(r(i,mm)/r(i,j))+cd*beta(i,j));
    end
    c1=sin(alpha)*sth(i)+cos(alpha)*cth(i);
    Vt(i)=V_inf*c1+sum1+sum2;
end

end

function z = simps(x,y,dim)
% got this from the mathworks website. credits below

%SIMPS  Simpson's numerical integration.
%   The Simpson's rule for integration uses parabolic arcs instead of the
%   straight lines used in the trapezoidal rule.
%
%   Z = SIMPS(Y) computes an approximation of the integral of Y via the
%   Simpson's method (with unit spacing). To compute the integral for
%   spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, SIMPS(Y) is the integral of Y. For matrices, SIMPS(Y) is a
%   row vector with the integral over each column. For N-D arrays, SIMPS(Y)
%   works across the first non-singleton dimension.
%
%   Z = SIMPS(X,Y) computes the integral of Y with respect to X using the
%   Simpson's rule. X and Y must be vectors of the same length, or X must
%   be a column vector and Y an array whose first non-singleton dimension
%   is length(X). SIMPS operates along this dimension.
%
%   Z = SIMPS(X,Y,DIM) or SIMPS(Y,DIM) integrates across dimension DIM of
%   Y. The length of X must be the same as size(Y,DIM).
%
%   Examples:
%   --------
%   % The integration of sin(x) on [0,pi] is 2
%   % Let us compare TRAPZ and SIMPS
%   x = linspace(0,pi,6);
%   y = sin(x);
%   trapz(x,y) % returns 1.9338
%   simps(x,y) % returns 2.0071
%
%   If Y = [0 1 2
%           3 4 5
%           6 7 8]
%   then simps(Y,1) is [6 8 10] and simps(Y,2) is [2; 8; 14]
%
%   -- Damien Garcia -- 08/2007, revised 11/2009
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>
%
%   See also CUMSIMPS, TRAPZ, QUAD.

%   Adapted from TRAPZ

%--   Make sure x and y are column vectors, or y is a matrix.
perm = []; nshifts = 0;
if nargin == 3 % simps(x,y,dim)
    perm = [dim:max(ndims(y),dim) 1:dim-1];
    yp = permute(y,perm);
    [m,n] = size(yp);
elseif nargin==2 && isscalar(y) % simps(y,dim)
    dim = y; y = x;
    perm = [dim:max(ndims(y),dim) 1:dim-1];
    yp = permute(y,perm);
    [m,n] = size(yp);
    x = 1:m;
else % simps(y) or simps(x,y)
    if nargin < 2, y = x; end
    [yp,nshifts] = shiftdim(y);
    [m,n] = size(yp);
    if nargin < 2, x = 1:m; end
end
x = x(:);
if length(x) ~= m
    if isempty(perm) % dim argument not given
        error('MATLAB:simps:LengthXmismatchY',...
            'LENGTH(X) must equal the length of the first non-singleton dimension of Y.');
    else
        error('MATLAB:simps:LengthXmismatchY',...
            'LENGTH(X) must equal the length of the DIM''th dimension of Y.');
    end
end

%-- The output size for [] is a special case when DIM is not given.
if isempty(perm) && isequal(y,[])
    z = zeros(1,class(y));
    return
end

%-- Use TRAPZ if m<3
if m<3
    if exist('dim','var')
        z = trapz(x,y,dim);
    else
        z = trapz(x,y);
    end
    return
end

%-- Simpson's rule
y = yp;
clear yp

dx = repmat(diff(x,1,1),1,n);
dx1 = dx(1:end-1,:);
dx2 = dx(2:end,:);

alpha = (dx1+dx2)./dx1/6;
a0 = alpha.*(2*dx1-dx2);
a1 = alpha.*(dx1+dx2).^2./dx2;
a2 = alpha.*dx1./dx2.*(2*dx2-dx1);

z = sum(a0(1:2:end,:).*y(1:2:m-2,:) +...
    a1(1:2:end,:).*y(2:2:m-1,:) +...
    a2(1:2:end,:).*y(3:2:m,:),1);

if rem(m,2) == 0 % Adjusting if length(x) is even
    state0 = warning('query','MATLAB:nearlySingularMatrix');
    state0 = state0.state;
    warning('off','MATLAB:nearlySingularMatrix')
    C = vander(x(end-2:end))\y(end-2:end,:);
    z = z + C(1,:).*(x(end,:).^3-x(end-1,:).^3)/3 +...
        C(2,:).*(x(end,:).^2-x(end-1,:).^2)/2 +...
        C(3,:).*dx(end,:);
    warning(state0,'MATLAB:nearlySingularMatrix')
end

%-- Resizing
siz = size(y); siz(1) = 1;
z = reshape(z,[ones(1,nshifts),siz]);
if ~isempty(perm), z = ipermute(z,perm); end

end
