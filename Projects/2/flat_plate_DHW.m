clear; close all; clc

%% location and radiation data

phi = 46*pi/180; % latitude (Trento)
n =  167; %15; % 15 June
delta = 23.45*pi/180*sin(2*pi*(n+284)/365); % declination 
omega_s = acos(-tan(delta)*tan(phi)); % sunset angle
t_sr = -12/pi*omega_s + 12; % sunrise time
t_ss = 12/pi*omega_s + 12; % sunset time
t =linspace(t_sr,t_ss,200)'; % time vector
t = t(2:end-1);
omega = pi/12*(t-12); % vector of hour angles 

Hh_day =  6.47e3; %1.65e3; % % Wh/m^2 daily global irradiation on a horizontal surface  

% ambient temperature
Ta_min = 16; %4; % °C 
Ta_max = 23; %8; % °C 
Ta = (Ta_min+Ta_max)/2 + (-Ta_min+Ta_max)/2*cos(2*pi/24*(t-12));   % temperature profile 

figure(1)
hold on
plot(t,Ta,'k')
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('T_a (°C)','fontsize',12,'fontname','times new roman')
set(gca,'fontsize',12,'fontname','times new roman')
set(gcf,'color','w')

%% flat plate collector features

A = 1.5; % m^2
eta0 = 0.8;
U = 5; % W/(m^2*K)

% fixed in/out water temp
Tin = 35; % °C
Tout = 45; % °C
% mains'water temperature
Tmains = 10; % °C
c = 4200; % J/(kg*K) water specific heat capacity

figure(2)
hold on
fplot(@(x)(eta0-U*x),'k')
grid on
xlim([0 0.2])
ylim([0 0.9])
xlabel('\it(T_{in}-T_{amb})/I','fontsize',12,'fontname','times new roman')
ylabel('\it\eta','fontsize',12,'fontname','times new roman')
set(gca,'fontsize',12,'fontname','times new roman')
set(gcf,'color','w')

%% profile of hourly radiation - simplified method! 
% we assume that the hourly profile of the radiation follows the same trend
% as the extraterrestrial radiaiton (in fact, this is not exactly the case because
% of beam/diffuse contributions)

cos_thetaz = sin(delta)*sin(phi) + cos(delta)*cos(phi)*cos(omega);
Iavg = Hh_day./trapz(t,(cos_thetaz)); % W/m^2
I = Iavg * (cos_thetaz); 
% trapz(t,I)

figure(3)
hold on
plot(t,I,'k')
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('I (W/m^2)','fontsize',12,'fontname','times new roman')
set(gca,'fontsize',12,'fontname','times new roman')
set(gcf,'color','w')
ylim([0 800])

%% Collector productivity

% collector hourly efficiency
eta = eta0 - U*(Tin-Ta)./I;

figure(4)
hold on
plot(t,eta,'k')
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('\it\eta','fontsize',12,'fontname','times new roman')
set(gca,'fontsize',12,'fontname','times new roman')
set(gcf,'color','w')
ylim([0 0.8])

ind = find(eta>0);

%  profile of harvested power
figure(5)
hold on
plot(t(ind),I(ind).*eta(ind)*A,'k')
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('Power (W)','fontsize',12,'fontname','times new roman')
set(gca,'fontsize',12,'fontname','times new roman')
set(gcf,'color','w')

% energy positively harvested (cumulative)
Energy = 3600*cumtrapz(t(ind),I(ind).*eta(ind)*A); % J 

% Amount of water heated up (from Tin to Tout)
Mh2o = Energy(end)/(c*(Tout-Tin))

% Amount of water heated up (from Tmains to Tout)
Mh2o_m = Energy(end)/(c*(Tout-Tmains))
