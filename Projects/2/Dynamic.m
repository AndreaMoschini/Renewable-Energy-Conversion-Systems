clear;
clc;
close all;






% Consumption profiles


% SELECT SEASON:
%     1 = winter (15 january)
%     2 = summmer (15 june)

season = 1;

if (season ~= 1 || season ~= 2)
    season = 1;
end


% we will assume that maximum temperature during the day is reached at t=tmax 
tmax = 14*60;


time = (0:60*24)'; % min
flow_rate_consumption = 0*time; % l/min
temperature_req = 0*time; % °C


if season == 1 % winter

    Tmains = 12; % ° C

    flow_rate_consumption(60*7:60*7+8) = 10; % l/min shower 1
    temperature_req(60*7:60*7+8) = 38; % °C

    flow_rate_consumption(60*7.25:60*7.25+2) = 6; % l/min face washing 1
    temperature_req(60*7.25:60*7.25+2) = 35; % °C

    flow_rate_consumption(60*7.5:60*7.5+3) = 6; % l/min breakfast wash
    temperature_req(60*7.5:60*7.5+3) = 32; % °C

    flow_rate_consumption(60*13:60*13+3) = 6; % l/min lunch wash
    temperature_req(60*13:60*13+3) = 40; % °C

    flow_rate_consumption(60*16:60*16+1) = 6; % l/min face washing 2
    temperature_req(60*16:60*16+1) = 35; % °C

    flow_rate_consumption(60*20:60*20+5) = 6; % l/min dinner wash
    temperature_req(60*20:60*20+5) = 42; % l/min dinner wash

    flow_rate_consumption(60*21.5:60*21.5+7) = 10; % l/min shower 2
    temperature_req(60*21.5:60*21.5+7) = 38; % °C

    n = 15;
    H = 3.6*1.649;     % [MJ/m^2]  daily global irradiation
    D = 3.6*0.713;     % [MJ/m^2]  daily diffuse irradiation 
    Tmin = 4;          % [°C] average monthly minimum temperature
    Tmax = 8;          % [°C] average monthly maximum temperature


else % summer

    Tmains = 16; 

    flow_rate_consumption(60*6:60*6+8) = 10; % l/min shower 1
    temperature_req(60*6:60*6+8) = 36; % ° C

    flow_rate_consumption(60*6.25:60*6.25+4) = 8; % l/min face washing 1
    temperature_req(60*6.25:60*6.25+4) = 32; % ° C

    flow_rate_consumption(60*6.5:60*6.5+3) = 6; % l/min breakfast wash
    temperature_req(60*6.5:60*6.5+3) = 32; % ° C

    flow_rate_consumption(60*12:60*12+3) = 6; % l/min lunch wash
    temperature_req(60*12:60*12+3) = 40; % ° C

    flow_rate_consumption(60*18:60*18+10) = 10; % l/min evenign shower
    temperature_req(60*18:60*18+10) = 35; % ° C

    flow_rate_consumption(60*19:60*19+5) = 6; % l/min dinner wash
    temperature_req(60*19:60*19+5) = 40; % l/min dinner wash

    flow_rate_consumption(60*20.5:60*20.5+8) = 10; % l/min shower 2
    temperature_req(60*20.5:60*20.5+8) = 36; % ° C

    n = 166;
    H = 3.6*6.471;     % [MJ/m^2]  daily global irradiation
    D = 3.6*2.528;     % [MJ/m^2]  daily diffuse irradiation 
    Tmin = 16;          % [°C] average monthly minimum temperature
    Tmax = 23;          % [°C] average monthly maximum temperature

end

figure(1)
ha(1)=subplot(2,1,1);
hold on
grid on
plot(time/60,flow_rate_consumption,'k')
ylabel('Flow rate (l/min)')

ha(2)=subplot(2,1,2);
hold on
grid on
plot(time/60,temperature_req,'k')
ylabel('Temp (°C)')
xlabel('hour')

linkaxes(ha,'x')



% METEOROLOGICAL INPUT 
lat=46.07; %latitude in deg
az=0; % panel azimuth in deg (facing south)
tilt=30; % tilt angle in deg
rho=0.2; % substrate reflectance


I0 = 1367; %[W/m^2] solar constant
phi = pi/180*lat; % latit in radianti
gamma = az*pi/180; % azimuth in rad % ELIMINA 
beta = pi/180*tilt; % tilt in rad
DT = Tmax-Tmin; % daily temperature variation

delta = 23.45*pi/180*sin(2*pi*(284+n)/365.25); % declination [rad]
omega_s = acos(-tan(phi)*tan(delta)); % sunrise/sunset angle on ground [rad]
t_sr = round((-12/pi*omega_s + 12)*60); % sunset time (m) rounded
t_ss = round((12/pi*omega_s + 12)*60); % sunrise time (m) rounded


N=1;             % number of modules installed
A = 2.14;        %[m^2] effective absorber surface
V = 1.7;         % [L] collector liquid content
alpha = 0.95;    % collector absorption

eta_zero = 0.802;   %optical efficiency at absorber
U = 4.3;            % [W / m^2K]
a1 = 4.28;          % [W / m^2K] coefficients of heat loss
a2 = 0.0064;        % [W / m^2K^2] coefficients of heat loss
flow_rate = 64.2;   % [L/h] recommended flow rate
epsilon = 0.6;   % heat exchanger effectiveness
tank_loss_coeff = 2.5; % [W / m^2K] tank loss coefficient
Thome = 18;         % [°C] room remperature
Tank_volume = 200;     % [kg] tank content mass

c_joule = 4200; % J/(kg*K) water specific heat capacity
c = c_joule/60; % W*min/(kg*K) water specific heat capacity
T_in = 30; % °C
water_density = 1;     % [kg/L] water density, assumed constant
m_dot = (flow_rate/60)*water_density;    % [kg/min] mass flow rate










% IRRADIATION CALCULATION

a_=0.409+0.5016*sin(omega_s-pi/3);
b_=0.6609-0.4767*sin(omega_s-pi/3);




% Ht = zeros(n_samples, 12);
% Dt = zeros(n_samples, 12);
% Bt = zeros(n_samples, 12);
% Hbt = zeros(n_samples, 12);
% Ta = zeros(n_samples, 12);
% solar_power = zeros(n_samples, 12);
% collector_power = zeros(n_samples, 12);
% deltaT = zeros(n_samples, 12);  %Tout - Tin (temperature)
% time = zeros(n_samples, 12);
% Tout = zeros(n_samples, 12);



% create timeseries that goes from sunrise to sunset
time_sun = (t_sr:t_ss)';
omega = pi/(12*60)*(time_sun-(12*60)); % hour angle of the sun
dt = time_sun(2)-time_sun(1);
domega = omega(2)-omega(1);


% ambient temperature profile
Ta = (Tmax+Tmin)/2 +(Tmax-Tmin)/2*cos(2*pi/(24*60)*(time-tmax)); % [°C] hourly ambient temperature


figure(2)
plot(time, Ta)
xlabel('Time [m]')
ylabel('Temperature [°C]')
title('Daily ambient temperature')


% calculate the hourly irradiation values. The total irradiation for
% each hour is calculated by multiplying these ratios for the average
% day irradiation (same goes for diffused component), then we calculate
% the beam component as difference
rd = pi/24*(cos(omega)-cos(omega_s))/(sin(omega_s)-omega_s*cos(omega_s));
rt = rd.*(a_+b_*cos(omega));
Dt = rd*D*dt/1;    % [MJ/m^2] in a time interval dt
Ht = rt*H*dt/1;    % [MJ/m^2] in a time interval dt
Bt = Ht - Dt;         % [MJ/m^2] in a time interval dt


% Check if correct
D_check = sum(Dt); % should be the same as D
H_check = sum(Ht); % should be the same as H

%Now that we calculated the hourly irradiation, we need to compute how
%much of it goes on our tilted surface

% Compute: thetai (angle of incidence of the sun rays on tilted plane)
% and thetaz (solar zenith angle)

cosThetaz=cos(omega)*cos(delta)*cos(phi)+sin(delta)*sin(phi);
cosThetai = cos(beta)*(sin(delta)*sin(phi) + cos(delta)*cos(omega)*cos(phi)) - sin(beta)*(cos(gamma)*(cos(phi)*sin(delta) - cos(delta)*cos(omega)*sin(phi)) - cos(delta)*sin(gamma)*sin(omega));

Rt=cosThetai./cosThetaz;


Hbt = (Bt.*Rt+Dt.*((1+cos(beta))/2)+((1-cos(beta))/2)*rho.*Ht)*1e6/(3600*dt); % [W/m^2] radiation (instantaneous power, function fo time)
%Hbt(:, i) = (Bt(:, i).*Rt+Dt(:, i).*((1+cos(beta))/2)+((1-cos(beta))/2)*rho.*Ht(:, i))*1e6*dt/(3600); % [W/m^2] radiation (instantaneous power, function fo time)
%Hbt(1) = 0; Hbt(end) = 0;
Hbt(Hbt<0) = 0;
% Hbt_(cosThetai<0 | cosThetaz<0)= 0;

Irr = zeros(length(time), 1);
Irr(time_sun) = Hbt;




%Plots
figure(3)
plot(time_sun, Ht, 'Color','b', 'LineStyle','-')
hold on
%plot(time, Hbt, 'Color','b', 'LineStyle','--')
plot(time_sun, Dt, 'Color','r', 'LineStyle','-')
plot(time_sun, Bt, 'Color','g', 'LineStyle','-')
xlabel('Time [m]')
xlim([0 24*60])
ylabel('Irradiation [MJ/m^2]')
legend('Total irradiation (H)', 'Diffused irradiation (D)', 'Beam component (B)')
title('Daily irradiation components ')
hold off


figure(4)
plot(time, Irr)
xlabel('Time [m]')
ylabel('Irradiation [W/m^2]')
title('Incident radiation on panel')


% SIM setup

I_timeseries = timeseries(Irr, time);
Ta_timeseries = timeseries(Ta, time);

mu_timeseries = timeseries(flow_rate_consumption, time);
Tu_timeseries = timeseries(temperature_req, time);

Ts0 = 30;       % [C] initial tank temperature
% we assume that the input tempertature in the collector is the same as the
% tank temperature
Tin_0 = Ts0;     % [C] initial collector input temperature
E_solar_0 = 0;
E_load_0 = 0;
E_boiler_0 = 0;
E_request_0 = 0;

%out = sim("Model\System_model.slx");


%% SIM - n days

% i separated E_load from E_request because E_load not always considers the
% difference of temperature between the request and the mains.





n_days = 5;

time_long = (0:(length(time)*n_days -1))';
Ts = zeros(length(time),n_days);
Tin = zeros(length(time),n_days);
Tout = zeros(length(time),n_days);
E_boiler = zeros(length(time),n_days);
E_solar = zeros(length(time),n_days);
E_load = zeros(length(time),n_days);
E_request = zeros(length(time),n_days);
Q_boiler = zeros(length(time),n_days);
Q_solar = zeros(length(time),n_days);
Q_load = zeros(length(time),n_days);
collector_pump = zeros(length(time),n_days);



for i=1:n_days

    out = sim("Model\System_model.slx");

    Ts(:, i) = out.Ts;
    Tin(:, i) = out.Tin;
    Tout(:, i) = out.Tout;
    E_boiler(:, i) = out.E_boiler;
    E_solar(:, i) = out.E_solar;
    E_load(:, i) = out.E_load;
    E_request(:, i) = out.E_request;
    Q_boiler(:, i) = out.Q_boiler;
    Q_solar(:, i) = out.Q_solar;
    Q_load(:, i) = out.Q_load;


% put final conditions of previous simulation as initial conditions of next
% simulation

    Ts0 = out.Ts(end);
    Tin_0 = out.Tin(end);
    E_solar_0 = out.E_solar(end);
    E_load_0 = out.E_load(end);
    E_boiler_0 = out.E_boiler(end);
    E_request_0 = out.E_request(end);

end

Ts = Ts(:);
Tin = Tin(:);
Tout = Tout(:);
E_boiler = E_boiler(:);
E_solar = E_solar(:);
E_load = E_load(:);
E_request = E_request(:);
Q_boiler = Q_boiler(:);
Q_solar = Q_solar(:);
Q_load = Q_load(:);
collector_pump = collector_pump(:);

% Plots


figure(5)
plot(time_long, Ts, 'Color', 'g')
hold on
plot(time_long, Tin, 'Color', 'b')
plot(time_long, Tout, 'Color', 'r')
xlabel('Time [min]')
ylabel('Temperature [°C]')
title('Collector and Tank temperatures')
legend('Tank', 'Collector IN', 'Collector OUT')
hold off



figure(6)
plot(time_long, E_request, 'Color', 'g')
hold on
plot(time_long, E_load, 'Color', 'y')
plot(time_long, E_boiler, 'Color', 'b')
plot(time_long, E_solar, 'Color', 'r')
xlabel('Time [min]')
ylabel('Energy [W*min]')
title('Energy share comparison')
legend('Requested', 'Load', 'Boiler', 'Solar')
hold off

% if we look closely, E_request is almost the sum between E_boiler and
% E_solar, representing the total energy requested by the whole system.
% However if we start with the initial Ts at 30 C we see that there is
% still a little part of energy missing in E_request to reach the sum of
% E_boiler with E_solar.
% This can be explained by the fact that in the tank model we have some
% energy losses that however are not accounted for in the total requested
% energy, but the solar+boiler energy still have to cover.
% NB: solar+boiler energy is higher than request energy

% If we cancel the tank evergy loss term we see that the trend inverts and
% the boiler+solar energy is lower than the request energy.
% This might be due to the fact that we are starting with Ts=30 C, which is
% a bit higher than the value at which it settles after a few days.
% If we have as initial Tsa lower value (like 22), we can see how
% solar+boiler energy gets much closer to the requested energy.



figure(7)
plot(time_long, Q_load, 'Color', 'g')
hold on
plot(time_long, Q_boiler, 'Color', 'b')
plot(time_long, Q_solar, 'Color', 'r')
xlabel('Time [min]')
ylabel('Power [W]')
title('Power share comparison')
legend('Load', 'Boiler', 'Solar')
hold off


% let's consider 2 energy inputs to the system (boiler and sun) and
% calculate how much the solar energy covers this request
eff = E_solar(end)/(E_solar(end)+E_boiler(end));


disp(['The solar energy covers the ', num2str(eff*100),' % of the total requested energy in the season ', num2str(season)]);

% NB: to find W*h instead of W*min, just divide the value by 60






%% SIM - sensitivity plots



n_days = 5;

time_long = (0:(length(time)*n_days -1))';
Ts = zeros(length(time),n_days);
Tin = zeros(length(time),n_days);
Tout = zeros(length(time),n_days);
E_boiler = zeros(length(time),n_days);
E_solar = zeros(length(time),n_days);
E_load = zeros(length(time),n_days);
E_request = zeros(length(time),n_days);
Q_boiler = zeros(length(time),n_days);
Q_solar = zeros(length(time),n_days);
Q_load = zeros(length(time),n_days);
collector_pump = zeros(length(time),n_days);



for i=1:n_days

    simIn = Simulink.SimulationInput("Model\System_model.slx");
    simIn = simIn.setVariable("epsilon", (epsilon));
    sinIn = simIn.setVariable("Tank_volume", (Tank_volume+0.2*Tank_volume));
    out = sim(simIn);

    Ts(:, i) = out.Ts;
    Tin(:, i) = out.Tin;
    Tout(:, i) = out.Tout;
    E_boiler(:, i) = out.E_boiler;
    E_solar(:, i) = out.E_solar;
    E_load(:, i) = out.E_load;
    E_request(:, i) = out.E_request;
    Q_boiler(:, i) = out.Q_boiler;
    Q_solar(:, i) = out.Q_solar;
    Q_load(:, i) = out.Q_load;


% put final conditions of previous simulation as initial conditions of next
% simulation

    Ts0 = out.Ts(end);
    Tin_0 = out.Tin(end);
    E_solar_0 = out.E_solar(end);
    E_load_0 = out.E_load(end);
    E_boiler_0 = out.E_boiler(end);
    E_request_0 = out.E_request(end);

end

Ts = Ts(:);
Tin = Tin(:);
Tout = Tout(:);
E_boiler = E_boiler(:);
E_solar = E_solar(:);
E_load = E_load(:);
E_request = E_request(:);
Q_boiler = Q_boiler(:);
Q_solar = Q_solar(:);
Q_load = Q_load(:);
collector_pump = collector_pump(:);

% Plots


figure(5)
plot(time_long, Ts, 'Color', 'g')
hold on
plot(time_long, Tin, 'Color', 'b')
plot(time_long, Tout, 'Color', 'r')
xlabel('Time [min]')
ylabel('Temperature [°C]')
title('Collector and Tank temperatures')
legend('Tank', 'Collector IN', 'Collector OUT')
hold off



figure(6)
plot(time_long, E_request, 'Color', 'g')
hold on
plot(time_long, E_load, 'Color', 'y')
plot(time_long, E_boiler, 'Color', 'b')
plot(time_long, E_solar, 'Color', 'r')
xlabel('Time [min]')
ylabel('Energy [W*min]')
title('Energy share comparison')
legend('Requested', 'Load', 'Boiler', 'Solar')
hold off




figure(7)
plot(time_long, Q_load, 'Color', 'g')
hold on
plot(time_long, Q_boiler, 'Color', 'b')
plot(time_long, Q_solar, 'Color', 'r')
xlabel('Time [min]')
ylabel('Power [W]')
title('Power share comparison')
legend('Load', 'Boiler', 'Solar')
hold off











%% SIM - sensitivity analysis but automated

% The simulation is still done over a few days to observe the influence of
% some parameters

% Parameters to change:
epsilon = 0.6;   % heat exchanger effectiveness
Tank_volume = 200;     % [kg] tank content mass



% Reference

eff_ref = parametric_analysis(epsilon, Tank_volume, time);

% epsilon ref +- 10%

eff_epsi1 = parametric_analysis((epsilon+0.1*epsilon), Tank_volume, time);
eff_epsi2 = parametric_analysis((epsilon-0.1*epsilon), Tank_volume, time);

% Tank volume +- 50%

% eff_tank1 = parametric_analysis(epsilon, (Tank_volume+0.5*Tank_volume), time);
% eff_tank2 = parametric_analysis(epsilon, (Tank_volume-0.5*Tank_volume), time);

% since we are calculating the efficiency as rapporto between energies,
% tank size doesn't matter i guess?



% Temp min and max +- 10%

T_min = Tmin + 0.1*Tmin;
T_max = Tmax;
Ta = (T_max+T_min)/2 +(T_max-T_min)/2*cos(2*pi/(24*60)*(time-tmax)); % [°C] hourly ambient temperature
Ta_timeseries = timeseries(Ta, time);
eff_Tmin1 = parametric_analysis(epsilon, Tank_volume, time);


T_min = Tmin - 0.1*Tmin;
T_max = Tmax;
Ta = (T_max+T_min)/2 +(T_max-T_min)/2*cos(2*pi/(24*60)*(time-tmax)); % [°C] hourly ambient temperature
Ta_timeseries = timeseries(Ta, time);
eff_Tmin2 = parametric_analysis(epsilon, Tank_volume, time);


T_min = Tmin;
T_max = Tmax + 0.1*Tmax;
Ta = (T_max+T_min)/2 +(T_max-T_min)/2*cos(2*pi/(24*60)*(time-tmax)); % [°C] hourly ambient temperature
Ta_timeseries = timeseries(Ta, time);
eff_Tmax1 = parametric_analysis(epsilon, Tank_volume, time);


T_min = Tmin;
T_max = Tmax - 0.1*Tmax;
Ta = (T_max+T_min)/2 +(T_max-T_min)/2*cos(2*pi/(24*60)*(time-tmax)); % [°C] hourly ambient temperature
Ta_timeseries = timeseries(Ta, time);
eff_Tmax2 = parametric_analysis(epsilon, Tank_volume, time);


%not showing tank size cause it doesn't matter
names = ["epsilon -10%", "Tmax -10%", "Tmin -10%", "Ref", "Tmin +10%", "Tmax +10%", "epsilon +10%"];
values = [eff_epsi2, eff_Tmax2, eff_Tmin2, eff_ref, eff_Tmin1, eff_Tmax1, eff_epsi1];

%names = ["epsilon -10%", "Tmin -10%", "Tmax -10%", "Ref", "Tmax +10%", "Tmin +10%", "epsilon +10%"];
%values = [eff_epsi2, eff_Tmin2, eff_Tmax2, eff_ref, eff_Tmax1, eff_Tmin1, eff_epsi1];

figure(8)
bar(names, values)
xlabel('Parameter')
ylabel('Efficiency')
title('Sensitivity analysis')