% EECS725 Homework 4

clear;
close all;

% Universal constants
c         = 3e8;       % speed of light (m/s)

% Problem constraints
f_start   = 100e6;      % chirp start frequency (hz)
f_end     = 300e6;      % chirp end frequency (hz)
tau       = 10e-6;      % chirp duration (s)
A_t1      = 0.1;        % amplitude of return, target 1
A_t2      = 1.0;        % amplitude of return, target 2
A_t3      = 0.25;       % amplitude of return, target 3
d_t1      = 25;         % delay of return, target 1 (ns)
d_t2      = 150;        % delay of return, target 2 (ns)
d_t3      = 160;        % delay of return, target 3 (ns)

% Simulation parameters
t_samp    = 1e-9;       % simulation sampling period (s)

%% Start simulation

% Generate chirp burst
k = (f_end-f_start) / tau; % chirp rate (hz/s)

t_burst = linspace(0,tau,(tau/t_samp)-1);
s       = cos( 2*pi*(f_start.*t_burst + 0.5*k.*(t_burst.^2)) )';
[Pss,w] = periodogram(s,ones(8192,1),8192);
f_plot  = w / (2*pi*t_samp);

% Generate composite radar return signal
r = zeros(2*length(s),1);
r_t1 = r;
r_t2 = r;
r_t3 = r;
r_t1(d_t1:length(s)+d_t1-1) = A_t1 * s;
r_t2(d_t2:length(s)+d_t2-1) = A_t2 * s;
r_t3(d_t3:length(s)+d_t3-1) = A_t3 * s;
r = r_t1 + r_t2 + r_t3;

% Multiply transmit and receive signals
mixer_out = s .* r(1:length(s));
mixer_psd_out = zeros(length(s),1);
mixer_psd_out(5000:5000+length(mixer_out)-1) = mixer_out;

% Spectral analysis
[Pmm,w] = periodogram(mixer_psd_out); %,ones(8192,1),8192);
f_mixer_plot = w / (2*pi*t_samp);

%% Plots-----------
figure(1)
plot(f_plot,10*log10(Pss));
title('Chirp PSD');
xlabel('Frequency (hz)');
ylabel('Power (dB)');
axis([50e6 400e6 -40 0]);
grid on;


figure(2)
plot(f_mixer_plot,10*log10(Pmm));
title('Return PSD');
xlabel('Frequency (hz)');
ylabel('Power (dB)');
axis([25e3 10e6 -40 40]);
grid on;