function x = S_fun(sig_speed, freqq)
    
    x = (0.1^2 * sig_speed * 600) ./ ((1+ 1.5.*((freqq.*600)./sig_speed)).^(5/3));

end