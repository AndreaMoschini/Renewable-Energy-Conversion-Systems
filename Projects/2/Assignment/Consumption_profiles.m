%% Consumption profiles

% !! Select a seasons (1=winter; 2=summer)
season = 1;

time = (0:60*24)'; % min
flow_rate_consumption = 0*time; % l/min
temperature_req = 0*time; % °C

if season == 2 % winter

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

end

figure(2)
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