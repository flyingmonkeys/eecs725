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

T = 2*(R*1000) / c; % round-time echo time (s)

v = sqrt( mu / (Re+h) );   % orbital velocity (km/s)
v_g = (v * Re) / (Re + h); % ground velocity (km/s)

PRF_min = 2 * (v_g*1000) / l;

max_swath_width = (c/2) * ( (1/PRF_min) - (2*tau) );

