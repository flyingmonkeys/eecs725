% EECS725 Homework 2

clear;

% Universal constants
c  = 3e8;      % speed of light (m/s)
mu = 398600;   % standard gravitational parameter for Earth (km^3/s^2)
Re = 6378.145; % average Earth radius (km)

% Problem constraints
f       = 1e9; 		 % carrier frequency (hz)
l       = 5;   		 % antenna array length (m)
tau     = 50e-6; 	 % pulse duration (s)
theta_i = 37*pi/180; % mid-swath incidence angle (radians)
h       = 330;       % orbit altitude (km)
w_gr    = 50;		 % swath width (km)

w_r = w_gr * sin(theta_i); % slant-range swath width (km)
alpha_spread = w_gr / Re; % alpha_far - alpha_near (rad)

% Compute geometric measurements
gamma = asin( (Re/(Re+h)) * sin(theta_i) ); % (radians)
alpha = theta_i - gamma; % (radians)
fprintf(1,'Gamma = %0.2f deg, Alpha = %0.2f deg\n',gamma*180/pi,alpha*180/pi);

% Compute slant ranges
R  = sqrt( Re^2 + (Re+h)^2 - 2*Re*(Re+h)*cos(alpha) ); % slant range (km)
R1 = sqrt( Re^2 + (Re+h)^2 - 2*Re*(Re+h)*cos(alpha-(alpha_spread/2)) ); % near range (km)
R2 = sqrt( Re^2 + (Re+h)^2 - 2*Re*(Re+h)*cos(alpha+(alpha_spread/2)) ); % far range (km)
fprintf(1,'Slant range R = %0.2f km, R1 = %0.2f km, R2 = %0.2f km\n',R,R1,R2);

% Calculate echo times
T  = 2*(R*1000) / c;    % round-time echo time (s)
T1 = (2 * R1*1000) / c; % near-end echo (s)
T2 = (2 * R2*1000) / c; % far-end echo (s)

% Compute PRFmax based on the swath width constraint
delta_T = 2 * (w_r*1000) / c; % (s)
PRF_max_swath_limit = 1 / (2*tau + delta_T); % (hz)
fprintf(1,'Absolute PRFmax = %0.2f hz\n',PRF_max_swath_limit);

% Compute PRFmin based on Doppler requirements
v   = sqrt( mu / (Re+h) ); % orbital velocity (km/s)
v_g = (v * Re) / (Re + h); % ground velocity (km/s)
PRF_min_doppler = 2 * (v*1000) / l; % (hz)
fprintf(1,'PRFmin (from Doppler) = %0.2f hz\n',PRF_min_doppler);

% Compute PRF min/max ranges based on multiple pulses in-the-air
N = 2;
PRF_min_abs = 1 / (T1 - tau);
PRF_max     = 2 / (T2 + tau);
while( (N/(T2+tau)) < PRF_max_swath_limit ) % only iterate up to PRF_max
    PRF_min = (N-1) / (T1 - tau);
    PRF_max = N / (T2 + tau);
    fprintf(1,'N = %d:PRF_min = %0.2f hz PRF_max = %0.2f hz\n', N, PRF_min, PRF_max);
    N = N + 1;
end

% Problem 3 - Compute max swath width before no PRF ranges are allowed
w_r_max = (c/2) * ( (1/PRF_min) - (2*tau) ); % (m)
max_swath_width = w_r_max / sin(theta_i); % (m)
fprintf(1,'Max swath width = %0.2f km\n',max_swath_width/1000);

