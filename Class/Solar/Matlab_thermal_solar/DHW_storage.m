clear; close all; clc

%% thermal storage example 
% Problem 8.3.1 from Duffie & Beckmann, Solar Engineering of Thermal Processes

M = 1500; % kg rho*V storage capacity
c = 4200; % J/(kg*K) water specific heat capacity
hA = 11.1; % W/K h*A loss coefficient
Ts0 = 45; % °C initial tank temperature
Ta = 20; % °C ambient temperature

% energy supplied by collector (hour by hour)
Qd = zeros(12,1);
Qd(end-3:end) = 1e6*[21, 41, 60, 75]'; % J
dQd = Qd/3600; % W power

% loads (hour by hour)
Qr = 1e6*[12 12 11 11 13 14 18 21 20 20 18 16]';
dQr = Qr/3600; % W power

time = 0:3600:11*3600; % time vector in sec. (associated to loads/supplies)

%% dynamic simulation
sim('storage_dyn')

figure(1)
hold on
plot(Ts.Time/3600,Ts.data,'k')
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('Storage temperature, T_s (°C)','fontsize',12,'fontname','times new roman')
set(gca,'fontsize',12,'fontname','times new roman')
set(gcf,'color','w')
ylim([25 55])