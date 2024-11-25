%% WT main
clc; close all; clear all;

% Simulation time setup
T_sim = 600;              % sim time [s]
dT = 0.1;                     % fixed size time interval [s]
time = 0:dT:T_sim;          % time vector [s]

%% DTU 10 MW DATA

rho_air = 1.225;            % air density [kg/m3]
H = 119;                    % hub height [m]
R_rotor = 178.3/2;               % rotor radius [m]
A_rotor = pi*R_rotor^2;     % rotor area [m2]
V_rated = 11.4;               % rated wind speed [m/s]
V_ci = 4;                   % cut-in wind speed [m/s]
V_co = 25;
Jeq = 1.56e8;                % equivalent inertia [kg*m2]
b = 2e5;                  % damping (equivalent losses) coefficient [N*m*s/rad]
P_rated = 10e6;             % rated power [W]
pole_pairs = 320;           % poles pairs
magnets_flux = 19.49;       % magnets flux [Wb]
stat_resistance = 64e-3;    % stator resistance [ohm]

%% Optimal parameters

omega_opt = 0.776;



%% POWER PROFILE
% Draw the power curve turbine. Under the cut-in speed the
% power must be 0 as well as over the cut-off one. For Vw < V_rated the
% maximum available power depends on the optimal value of Cp 

V_w = (0:0.1:V_co+5)';


% working regions 
I1 = find(V_w<V_ci);
I2 = find(V_w>=V_ci & V_w<=V_co);



lambda_stall = (omega_opt*R_rotor)./V_w;
Cp_stall = PowerFactor(lambda_stall, 0);
Cp_stall(find(Cp_stall<0)) = 0;



Pw = 0*V_w;
Pw(I2) = 0.5*rho_air*V_w(I2).^3*A_rotor.*Cp_stall(I2);

figure;
hold on
plot(V_w,Pw/1e6,'r','linewidth',1.5)
grid on
xlabel('\itV_w\rm (m/s)','FontSize',12,'fontname','times new roman')
ylabel('\itP_r\rm (MW)','FontSize',12,'fontname','times new roman')
set(gca,'FontSize',12,'fontname','times new roman')

%% SIMULATION WITH RAMP WIND PROFILE 

% -- SETUP
V_end = V_co + 5; % final wind speed
Vw_ramp = linspace(1.1*V_ci,V_end,length(time));

Vw_timeseries = timeseries(Vw_ramp',time');
omega_rotor_0 = omega_opt; % initial condition (initial rotor speed)

figure(2)
hold on
plot(out.V_w,out.P_gen/1e6,'k')  %data coming from simulink simulation
plot(V_w,Pw/1e6,'r','linewidth',1.5)