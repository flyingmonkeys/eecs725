% EECS725 Homework 3

clear;
close all;

% Universal constants
c         = 3e8;       % speed of light (m/s)

% Problem constraints
h         = 3e3;       % height (m)
lambda    = 10e-2;     % radar wavelength (m)
beta      = 25*pi/180; % 3dB beam width of antenna (Gaussian beam shape)
Go_dB     = 16;        % antenna peak gain (dB)
tau       = 1e-9;      % pulse duration (s)
sigma_o_0 = 0.1;       % backscattering coefficient

Go        = 10^(Go_dB/10); % antenna gain (linear)

% Simulation parameters
t_max     = 1.3e-6;    % max simulation time
N         = 1000;      % simulation granularity
t         = linspace(0,t_max,N); % time

%% Start simulation

t0 = h/c; % time that front part of transmitted pulse hits the earth
t1 = t0 + (c*tau/2); % time that back part of transmitted pulse hits the earth

% Compute radii for leading and trailing edge of wave pulse
r2 = sqrt((h + c.*t./2).^2 - h^2); % (m)
r1 = sqrt((h + (c.*t./2) + (c*tau/2)).^2 - h^2); % (m)
r2(1) = 0; % trailing edge has not contacted ground yet, so first time unit is zero
R = sqrt(((r1+r2)/2).^2 + h^2);

% Illuminated area (annulus)
area = pi*(r1.^2 - r2.^2); % (m)

% Theta (angle of center of annulus with respect to nadir)
theta = asin( ((r1+r2)./2) ./ sqrt(((r1+r2)/2).^2 + h^2) ); % (rad)
    
% Antenna gain at theta
G = Go * exp(-2.773 * ((theta/beta).^2 + (theta/beta).^2)); % (linear)
    
% Terrain backscaterring coefficient at theta
sigma_o = sigma_o_0 * (cos(theta).^9); % (linear)
    
% Distance from annulus center back to transmitter
R = sqrt( h^2 + ((r1+r2)/2).^2 ); % (m)

% Pr/Pt
PrPt = (lambda^2 * G.^2 .* area .* sigma_o) ./ ((4*pi)^3 * R.^4); % (linear)

%% Plots---------------------

figure(1)
plot(t*1e6,area);
title('Illuminated Area vs. time');
xlabel('time (us)');
ylabel('Area (sq m)');
grid on;

figure(2)
plot(t*1e6,theta*180/pi);
title('Theta vs. time');
xlabel('time (us)');
ylabel('Theta (deg)');
grid on;

figure(3)
plot(t*1e6,10*log10(PrPt));
title('Pr/Pt vs. time');
xlabel('time (us)');
ylabel('Pr/Pt (dB)');
grid on;

figure(4)
plot(t*1e6,10*log10(G));
title('Gain vs. time');
grid on;
xlabel('time (us)');
ylabel('Gain (dB)');
% Contribution is squared, so multiply dB spread by 2

figure(5)
plot(t*1e6,10*log10(sigma_o));
title('Sigma vs. time');
xlabel('time (us)');
ylabel('sigma (dB)');
grid on;

figure(6)
plot(t*1e6,r1);
hold on;
plot(t*1e6,r2);
title('r1 and r2 vs. time');
xlabel('time (us)');
ylabel('Radius (m)');
grid on;