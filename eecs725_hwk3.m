% EECS725 Homework 3

clear;

% Universal constants
c  = 3e8;      % speed of light (m/s)

% Problem constraints
h         = 3e3;       % height (m)
lambda    = 10e-2;     % radar wavelength (m)
beta      = 25*pi/180; % 3dB beam width of antenna (Gaussian beam shape)
Go        = 16;        % antenna peak gain (dB)
tau       = 1e-9;      % pulse duration (s)
sigma_o_0 = 0.1;       % backscattering coefficient

t0 = h/c; % time that front part of transmitted pulse hits the earth
t1 = t0 + (c*tau/2); % time that back part of transmitted pulse hits the earth
%t = linspace(0,t0+1.3e-6,100); % time after the first transmitted energy hits the earth

idx = 1;
for time=0:12e-6/100:12e-6
    r1(idx) = sqrt(h^2 - (h - c*time)^2);
    r2(idx) = sqrt(h^2 - (h - c*(time-(c*tau/2)))^2);
    t(idx) = time;
    idx = idx + 1;
end

figure(1)
plot(t,r1);