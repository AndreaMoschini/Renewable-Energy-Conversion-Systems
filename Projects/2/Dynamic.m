clear;
clc;
close all;



%% Consumption profiles


% SELECT SEASON:
%     1 = winter (15 january)
%     2 = summmer (15 june)

season = 1;

if (season ~= 1 || season ~= 2)
    season = 1;
end


% we will assume that maximum temperature during the day is reached at t=tmax 
tmax = 14;


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



%% METEOROLOGICAL INPUT 
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
epsilon = 0.04;  % emissions
eta_zero = 0.802;   %optical efficiency at absorber
U = 4.3;            % [W / m^2K]
a1 = 4.28;          % [W / m^2K] coefficients of heat loss
a2 = 0.0064;        % [W / m^2K^2] coefficients of heat loss
flow_rate = 64.2;   % [L/h] recommended flow rate
effectiveness = 0.6;   % heat exchanger effectiveness
tank_loss_coeff = 2.5; % [W/K] tank loss coefficient
Thome = 18;         % [°C] room remperature

c_joule = 4200; % J/(kg*K) water specific heat capacity
c = c_joule/3600; % Wh/(kg*K) water specific heat capacity
T_in = 30; % °C
Tin = T_in + 273.15;   % K
water_density = 1;     % [kg/L] water density, assumed constant
m_dot = (flow_rate/60)*water_density;    % [kg/m] mass flow rate

%% IRRADIATION CALCULATION

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
Ta = (Tmax+Tmin)/2 +(Tmax-Tmin)/2*cos(2*pi/24*(time_sun-tmax)); % [°C] hourly ambient temperature


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
Hbt(1) = 0; Hbt(end) = 0;
Hbt(Hbt<0) = 0;
% Hbt_(cosThetai<0 | cosThetaz<0)= 0;



%Plots
figure(2)
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


%Now calculate

eta = eta_zero - U * (T_in - Ta) ./ Hbt;

ind = find(eta > 0);

%Total solar power on tilted panel

solar_power = Hbt.*A;

collector_power = eta.*Hbt.*A;
collector_power(collector_power<0) = 0;
collector_power(isnan(collector_power)) = 0;

%Knowing the power captured by the collector, the fluid flow rate (assumed constant at the recommended value) and
%the specific heat, we can calculate the temperature change of the
%fluid between the entrance and the exit of the collector. Keep in mind
%that the time is expressed in hours, therefore we multiply the flow
%rate by dt

deltaT = collector_power/(flow_rate*dt*c);
