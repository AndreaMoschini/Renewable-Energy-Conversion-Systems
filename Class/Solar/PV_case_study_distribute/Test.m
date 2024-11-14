clear;
close all;
clc

Tilt_v = linspace(0, 45, 50)';
AEP_tilt = 0*Tilt_v;

for kk = 1: length(Tilt_v)

    tilt = Tilt_v(kk);
    Fixed_PV_productivity
    close all; clc
end

plot(Tilt_v,AEP_tilt, 'k')
grid on