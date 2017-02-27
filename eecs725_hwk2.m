% EECS725 Homework 2

clear;

% Universal constants
c = 3e8;       % speed of light (m/s)
mu = 398600;   % standard gravitational parameter for Earth (km^3/s^2)
Re = 6378.145; % average Earth radius (km)

% Problem constraints
f       = 1e9; 		% carrier frequency (hz)
l       = 5;   		% antenna array length (m)
tau     = 50e-6; 	% pulse duration (s)
theta_i = 37*pi/180; 	% mid-swath incidence angle (radians)
h       = 330;          % orbit altitude (km)
w_gr    = 50;		% swath width (km)

w_r = w_gr * sin(theta_i); % slant-range swath width (km)

gamma = asin( (Re/(Re+h)) * sin(theta_i) ); % (radians)
alpha = theta_i - gamma; % (radians)

R = sqrt( Re^2 + (Re+h)^2 - 2*Re*(Re+h)*cos(alpha) ); % slant range (km)
R1 = R - (w_r/2); % near range (km)
R2 = R + (w_r/2); % far range (km)

T = 2*(R*1000) / c; % round-time echo time (s)
T1 = (2 * R1*1000) / c; % near-end echo (s)
T2 = (2 * R2*1000) / c; % far-end echo (s)

delta_T = 2 * (w_r*1000) / c; % (s)
PRF_max_swath_limit = 1 / (2*tau + delta_T); % (hz)
fprintf(1,'Absolute PRFmax = %0.2f hz\n',PRF_max_swath_limit);

v = sqrt( mu / (Re+h) );   % orbital velocity (km/s)
v_g = (v * Re) / (Re + h); % ground velocity (km/s)

%PRF_min = 2 * (v*1000) / l;
N = 2;
PRF_min_abs = 1 / (T1 - tau);
PRF_max = 2 / (T2 + tau);
while( (N/(T2+tau)) < PRF_max_swath_limit )
    PRF_min = (N-1) / (T1 - tau);
    PRF_max = N / (T2 + tau);
    fprintf(1,'PRF_min = %0.2f PRF_max = %0.2f\n', PRF_min, PRF_max);
    N = N + 1;
end

w_r_max = (c/2) * ( (1/PRF_min_abs) - (2*tau) ); % (m)
max_swath_width = w_r_max / sin(theta_i); % (m)
fprintf(1,'Max swath width = %0.2f km\n',max_swath_width/1000);

