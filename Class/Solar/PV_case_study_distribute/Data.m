n=[15
    45
    74
    105
    135
    166
    196
    227
    258
    288
    319
    349];  % ordinal number of XV day of each month

Month_day = [31
    28
    31
    30
    31
    30
    31
    31
    30
    31
    30
    31]; % duration (in days) of each month

% Radiation data for Trento 

H= 3.6*[1.649	
    2.514	
    3.889	
    4.903	
    5.636	
    6.471	
    6.488	
    5.575	
    4.313	
    2.812	
    1.653	
    1.277]; % [MJ/m^2]  daily global irradiation 

D = 3.6*[0.713	
    1.039	
    1.502	
    2.021	
    2.408	
    2.528	
    2.413	
    2.138	
    1.707	
    1.236	
    0.818	
    0.624]; % [MJ/m^2]  daily diffuse irradiation 

Tmin=[4  3  5 8 12 16 18  18 15  12 7 5]';  % [°C] average monthly minimum temperature
Tmax=[8   8  10 13    18  23  25  25  21  16   12 9]'; % [°C] average monthly maximum temperature

% we will assume that maximum temperature during the day is reached at t=tmax 
tmax = 14; 