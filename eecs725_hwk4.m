% EECS725 Homework 4

clear;
close all;

% Universal constants
c         = 3e8;       % speed of light (m/s)

% Problem constraints
f_start   = 100e6;      % chirp start frequency (hz)
f_end     = 300e6;      % chirp end frequency (hz)
tau       = 10e-6;      % chirp duration (s)

% Simulation parameters
t_samp    = 1e-9;       % simulation sampling period (s)

%% Start simulation

% Generate chirp burst
k = (f_end-f_start) / tau; % chirp rate (hz/s)

t_burst = linspace(0,tau,(tau/t_samp)-1);
s       = cos( 2*pi*(f_start.*t_burst + 0.5*k.*(t_burst.^2)) );
[Pss,w] = periodogram(s);
f_plot  = w / (2*pi*t_samp);

%% Plots-----------
figure(1)
plot(f_plot,10*log10(Pss));
title('Chirp PSD');
xlabel('Frequency (hz)');
ylabel('Power (dB)');
axis([50e6 400e6 -40 0]);
grid on;

