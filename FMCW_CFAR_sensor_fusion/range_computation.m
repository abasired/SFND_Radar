% Calculate range of a vehicle using beat frequency

f_beat = [0, 1.1e6, 13e6, 24e6];

% Speed of light
c = 3*10^8;

% Bandwidth of chirp for range resolution
delta_r = 1;
chirp_BW = c/(2*delta_r);

% Chirp time
range_max = 300;
Ts = 5.5 * (range_max*2/c);

% Range computation
calculated_range = c*Ts*f_beat/(2*chirp_BW);

disp(calculated_range);




