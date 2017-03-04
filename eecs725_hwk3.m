% EECS725 Homework 3

clear;
close all;

% Universal constants
c  = 3e8;      % speed of light (m/s)

% Problem constraints
h         = 3e3;       % height (m)
lambda    = 10e-2;     % radar wavelength (m)
beta      = 25*pi/180; % 3dB beam width of antenna (Gaussian beam shape)
Go        = 16;        % antenna peak gain (dB)
tau       = 1e-9;      % pulse duration (s)
sigma_o_0 = 0.1;       % backscattering coefficient

t_max = 1.3e-6;

t0 = h/c; % time that front part of transmitted pulse hits the earth
t1 = t0 + (c*tau/2); % time that back part of transmitted pulse hits the earth
%t = linspace(0,t0+1.3e-6,100); % time after the first transmitted energy hits the earth

N = 1000; % simulation granularity
r1 = zeros(N,1);
r2 = zeros(N,1);
t = zeros(N,1);
area = zeros(N,1);
theta = zeros(N,1);
PrPt = zeros(N,1);

idx = 1;
for time=0:t_max/N:t_max
    r1(idx) = sqrt(h^2 - (h - c*time)^2);
    if( time < tau )
        r2(idx) = 0;
    else
        r2(idx) = sqrt(h^2 - (h - (c*time) + (c*tau/2))^2);
    end
    t(idx) = time;
    
    % Area
    area(idx) = pi*(r1(idx)^2 - r2(idx)^2);
    
    % Theta
    theta(idx) = asin(((r1(idx)+r2(idx))/2)/h);
    
    % Antenna gain
    G = Go * exp(-2.773 * ((theta(idx)/beta)^2 + (theta(idx)/beta)^2));
    
    % Terrain backscaterring coefficient
    sigma_o = sigma_o_0 * (cos(theta(idx))^9);
    
    % Pr/Pt
    PrPt(idx) = (lambda^2 * G^2 * area(idx)) / ((4*pi)^3 * h);
    
    idx = idx + 1;
end

figure(1)
plot(t,r1);
hold on;
plot(t,r2);
title('r1 and r2 vs. time');

figure(2)
plot(t,area);
title('Area');

figure(4)
plot(t,theta*180/pi);
title('Theta vs. time');

figure(5)
plot(t,20*log10(PrPt));
title('Pr/Pt vs. time');