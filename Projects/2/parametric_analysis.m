function [efficiency] = parametric_analysis(Effectiveness,TankVol, time)


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
    simIn = simIn.setVariable("epsilon", Effectiveness);
    sinIn = simIn.setVariable("Tank_volume", TankVol);
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
efficiency = E_solar(end)/(E_solar(end)+E_boiler(end));



end

