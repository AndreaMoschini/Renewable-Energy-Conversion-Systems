clear; close all; clc
% Montly and Yearly productivity of a fixed PV module 

% METEOROLOGICAL INPUT 
lat=46.07; %latitude in deg
az=0; % panel azimuth in deg
tilt=20; % tilt angle in deg
rho=0.; % substrate reflectance
Data

% PV panel technical data 
N=1; % number of modules installed
Vmpp0=18.4; %[V] mpp voltage (@ 1000 W/m2)
Impp0=5.43; %[A] mpp current (@ 1000 W/m2)
Voc0=22.95; % [V] OC voltage (@ 1000 W/m2)
Isc0=5.85; % [A] SC current (@ 1000 W/m2)
chi_oc=-0.36; % [%/°C] Voc temperature coefficient
chi_sc = 0.02; % [%/°C] Isc temperature coefficient
NOCT=47; % [°C] normal operating cell temperature

I0=1367; %[W/m^2] solar constant
phi=pi/180*lat; % latit in radianti
gamma=az*pi/180; % azimuth in rad % ELIMINA 
beta=pi/180*tilt; % tilt in rad
DT=Tmax-Tmin; % daily temperature variation 
FF=(Vmpp0*Impp0)/(Voc0*Isc0);

delta=23.45*pi/180*sin(2*pi*(284+n)/365.25); % declination [rad]
omega_s=acos(-tan(phi)*tan(delta)); % sunrise/sunset angle on ground [rad]
t_sr = -12/pi*omega_s + 12; % sunset time (h) 
t_ss = 12/pi*omega_s + 12; % sunrise time (h) 

% calculation of diffuse component from global (cross-check) 
% H0=86400/(pi*10^6)*I0*(1+0.034*cos(2*pi*n/365.25)).*(cos(phi)*cos(delta).*sin(omega_s)+omega_s*sin(phi).*sin(delta)); % [MJ/m^2] extraterrestrial irradiance (monthly avg) 
% K=H./H0;
% D_k=H.*(1.39-4.027*K+5.553*K.^2-3.108*K.^3); %[MJ/m^2] estimated diffuse radiation component

% figure(1)
% hold on
% plot(1:12,D,'k')
% plot(1:12,D_k,'--r')
% grid on
% xlabel('time (h)','fontsize',12,'fontname','times new roman')
% ylabel('\itD\rm (MJ/m^2)','fontsize',12,'fontname','times new roman')
% legend('data','correlation')
% set(gca,'fontsize',12,'fontname','times new roman')
% set(gcf,'color','w')

% Calculation Of Hourly Irradiation Values
a_=0.409+0.5016*sin(omega_s-pi/3);
b_=0.6609-0.4767*sin(omega_s-pi/3);

% H_check = zeros(12,1);
% D_check = zeros(12,1);
Energy = zeros(12,1); 

for i=1:12 % loop over months
    time= linspace(t_sr(i),t_ss(i),500)';    
    omega = pi/12*(time-12); % hour angle of the sun 
    dt = time(2)-time(1);
    domega = omega(2)-omega(1);

    Ta = (Tmax(i)+Tmin(i))/2 +(Tmax(i)-Tmin(i))/2*cos(2*pi/24*(time-tmax)); % [°C] hourly ambient temperature

    rd = pi/24*(cos(omega)-cos(omega_s(i)))/(sin(omega_s(i))-omega_s(i)*cos(omega_s(i)));
    rt = rd.*(a_(i)+b_(i)*cos(omega));
    Dt = rd*D(i)*dt/1; % [MJ/m^2] in a time interval dt
    Ht = rt*H(i)*dt/1; % [MJ/m^2] in a time interval dt
    Bt = Ht - Dt; % [MJ/m^2] in a time interval dt

    % D_check(i) = sum(Dt); % should be the same as D
    % H_check(i) = sum(Ht); % should be the same as H
    
    cosThetaz=cos(omega)*cos(delta(i))*cos(phi)+sin(delta(i))*sin(phi);
    cosThetai = cos(beta)*(sin(delta(i))*sin(phi) + cos(delta(i))*cos(omega)*cos(phi)) - sin(beta)*(cos(gamma)*(cos(phi)*sin(delta(i)) - cos(delta(i))*cos(omega)*sin(phi)) - cos(delta(i))*sin(gamma)*sin(omega));

    Rt=cosThetai./cosThetaz;
    Hbt=(Bt.*Rt+Dt.*((1+cos(beta))/2)+((1-cos(beta))/2)*rho.*Ht)*1e6/(3600*dt); % [W/m^2] radiation (instantaneous power, function fo time)
    Hbt(1) = 0; Hbt(end) = 0; 
    Hbt(Hbt<0) = 0;
    % Hbt_(cosThetai<0 | cosThetaz<0)= 0;

    figure(2)
    hold on 
    plot(time,Hbt)

    Tm=Ta+(NOCT-20)*Hbt/800; % [°C] module temperature 
    Isc=Isc0*(1+chi_sc/100*(Tm-25)).*Hbt/1000; %[A] SC current
    Voc=Voc0*(1+chi_oc/100*(Tm-25)); % [V] OC voltage
    Power=N*Isc.*Voc*FF; % Watt 
    Energy(i) = trapz(time,Power); % Wh generated in the avg day of the current month

    figure(3)
    hold on 
    plot(time,Power)

end
figure(2)
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('Incident radiation (W/m^2)','fontsize',12,'fontname','times new roman')
legend('data','correlation')
set(gca,'fontsize',12,'fontname','times new roman')
legend('jan.','feb.','mar.','apr.','may','jun.','jul.','aug.','sep.','oct.','nov.','dic.')
ylim([0 1e3])

figure(3)
grid on
xlabel('time (h)','fontsize',12,'fontname','times new roman')
ylabel('Power output (W)','fontsize',12,'fontname','times new roman')
legend('data','correlation')
set(gca,'fontsize',12,'fontname','times new roman')
legend('jan.','feb.','mar.','apr.','may','jun.','jul.','aug.','sep.','oct.','nov.','dic.')
ylim([0 100])

AEP = Month_day'*Energy % annual energy production (Wh)
CF = AEP/(N*Vmpp0*Impp0*8760)