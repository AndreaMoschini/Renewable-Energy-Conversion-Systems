function Cp = PowerFactor(lambda,pitch_angle)

    % EXPONENTIAL Power factor
    one_over_kappa = 1./(lambda + 0.08*pitch_angle) - 0.035/(1+pitch_angle^3);
    Cp = 0.55 * (116 * one_over_kappa - 0.4*pitch_angle - 5) .* exp(-21*one_over_kappa);

    % SINUSOIDAL Power factor
    % Cp = (0.44 - 0.0167*pitch_angle) * sin( pi/2*(lambda-3)/(7.5-0.15*pitch_angle) ) - (lambda-3)*0.00184*pitch_angle;
end
