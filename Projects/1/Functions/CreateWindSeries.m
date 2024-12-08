function [v] = CreateWindSeries(time_series, Tsim, speed, frequency, phase_shift)

    v = zeros(length(time_series), 1);

    for i=1:1:length(time_series)

        v(i) = speed + sum(    (sqrt(  (2 .* KaimalSpectrum(speed, frequency))  ./  Tsim   )     .*  cos(2 .* pi .* frequency .* time_series(i) - phase_shift))    );

    end


end