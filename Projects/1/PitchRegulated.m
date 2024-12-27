%% WT main
clc; close all; clear all;

addpath("Functions\")
addpath("Models\")

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


%% POWER FACTOR MAP
% The map is used only for visualizing the plot, in the model the fitting
% formula is used for each time instant (NOT a look-up method)

% -- BUILD AXIS
theta_vec = (0:5:30)';                   % blades pitch angle axis [deg]
lambda_sweep = (0.1:0.1:18)';              % tip speed ratio axis []    

% -- CREATE MAP
Cp_map = zeros(length(lambda_sweep),length(theta_vec));
figure
hold on
for i=1:1:length(theta_vec)
    theta = theta_vec(i);
    % -- EXPONENTIAL FITTING FORMULA
    Cp_map(:,i) = PowerFactor(lambda_sweep,theta);

    % -- SINUSOIDAL FITTING FORMULA
    % Cp_map(:,i) = (0.44 - 0.0167*theta) * sin(pi/2*(lambda_sweep-3)/(7.5-0.15*theta))-(lambda_sweep-3)*0.00184*theta;
    plot(lambda_sweep,Cp_map(:,i),'-')
    leg{i} = ['\beta = ',num2str(theta),'°'];
end
% clear('one_over_beta');

% -- FIND PEAK POWER FACTOR VALUE lambda and pitch
lambda_opt = lambda_sweep(find(~(Cp_map(:,1) - max(Cp_map(:,1)))));     % lambda corresponding to the peak (pitch = 0°)
% one_over_beta_opt = 1./(lambda_opt) - 0.035;
% Cp_opt = 0.5 * (116 * one_over_beta_opt - 5)*exp(-21*one_over_beta_opt);

Cp_opt = max(Cp_map(:,1)); %(0.44) * sin(pi/2*(lambda_opt-3)/(7.5));                       % optimal Cp
% clear('one_over_beta_opt');
Trated = Cp_opt/lambda_opt*0.5*rho_air*A_rotor.*V_rated^2*R_rotor;

plot(lambda_opt,Cp_opt,'or')
hold off
ylim([0 0.5])
title('Power factor [-]')
xlabel(' \lambda [-]')
ylabel('C_p')
legend(leg)




%% POWER PROFILE
% Draw the power curve turbine. Under the cut-in speed the
% power must be 0 as well as over the cut-off one. For Vw < V_rated the
% maximum available power depends on the optimal value of Cp 

V_w = (0:0.1:V_co+5)';
Prated = 0.5*rho_air*V_rated^3*A_rotor*Cp_opt;

% working regions 
I1 = find(V_w<V_ci);
I2 = find(V_w>=V_ci & V_w<=V_rated);
I3 = find(V_w>V_rated & V_w<=V_co);
I4 = find(V_w>V_co);

Pw = 0*V_w;
Pw(I2) = 0.5*rho_air*V_w(I2).^3*A_rotor*Cp_opt;
Pw(I3) = Prated; %*ones(size(I3));

figure;
hold on
plot(V_w,Pw/1e6,'r','linewidth',1.5)
grid on
xlabel('\itV_w\rm (m/s)','FontSize',12,'fontname','times new roman')
ylabel('\itP_r\rm (MW)','FontSize',12,'fontname','times new roman')
set(gca,'FontSize',12,'fontname','times new roman')

%% SIMULATION WITH RAMP WIND PROFILE 

% -- SETUP
V_end = V_co +5; % final wind speed
Vw_ramp = linspace(1.1*V_ci,V_end,length(time));

Vw_timeseries = timeseries(Vw_ramp',time');
omega_rotor_0 = Vw_ramp(1)*lambda_opt/R_rotor; % initial condition (initial rotor speed)
Noise_switch = 0;

% Noise data, needed otherwise simulink throws an error. In this example
% it is not used anyway

noise_data = repelem(0, length(time));
noise_timeseries = timeseries(noise_data', time');

figure(2)
hold on
plot(out.V_w,out.P_gen/1e6,'k')  %data coming from simulink simulation
plot(V_w,Pw/1e6,'r','linewidth',1.5)