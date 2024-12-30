clear; close all; clc

% !! Remember to comment out lines 1 (clear, ...) and 7 (tilt = ...) in
% Fixed_PV_productivity.m before running this script !!

Tilt_v = linspace(0,45,50)'; 
AEP_tilt = 0*Tilt_v; 

for kk = 1: length(Tilt_v)
    
    tilt = Tilt_v(kk);
    Fixed_PV_productivity
    close all; clc

    AEP_tilt(kk) = AEP; 
end

plot(Tilt_v,AEP_tilt,'k')
grid on
xlabel('\beta (deg)')
ylabel('AEP (Wh)')

[AEP_max,kk_max] = max(AEP_tilt);
tilt_max = Tilt_v(kk_max)
