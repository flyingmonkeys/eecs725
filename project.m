% EECS725 Project

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
rho       = 64;         % pulse compression ratio

%% Start simulation

% Derived parameters
bw        = f_end - f_start; % bandwidth (hz)
t_b       = 1 / bw;          % chip period (s)
fc        = (f_start + f_end) / 2; % center (carrier) frequency for P4 (hz)

% Simulation parameters
f_samp    = 10;               % samples/chip
t_samp    = t_b / f_samp;    % simulation sampling period (s)

% Generate chirp burst
k = (f_end-f_start) / tau; % chirp rate (hz/s)

t_burst = linspace(0,tau,(tau/t_samp)-1);
s       = cos( 2*pi*(f_start.*t_burst + 0.5*k.*(t_burst.^2)) )';

% Get PSD of chirp
[Pss,w] = periodogram(s); %,ones(8192,1),8192);
f_plot  = w / (2*pi*t_samp);

% Get autocorrelation of chirp
%xc_lfm = xcorr2(s,9000);
xc_lfm = xcorr(s);

% Generate P4 signal (baseband phases)
theta = zeros(rho,1);
theta_abs(1) = 0;
for i=2:rho+1
    theta_abs(i) = pi * ( (((i-1)^2)/rho) - (i-1) );
theta(i-1) = theta_abs(i-1); % - theta_abs(i-1); 
%    theta(i) = 2*pi*i/rho;
end

% Generate P4 signal (passband signal)
s_p4 = zeros(length(t_burst),1);
pn = sign(rand-0.5);
for i=1:length(t_burst)
ss(i) = mod(floor(i/f_samp),rho)+1;
if(i > 1)
    if(ss(i) ~= ss(i-1))
        pn = sign(rand-0.5);
    end
end
%    s_p4_lp(i) = exp( j * theta(mod(floor(i/f_samp),rho)+1) );
%s_p4_lp(i) = 1; % DC test signal
s_p4_lp(i) = exp( j * ((pn+1)/2) );
%    s_p4(i) = real( s_p4_lp(i) * exp( j*2*pi*fc*(i-1)*t_samp ) );
%s_p4(i) = cos( 2*pi*fc*(i-1)*t_samp + (pi*((pn+1)/2)) );
s_p4(i) = cos( 2*pi*fc*(i-1)*t_samp + theta(mod(floor(i/f_samp),rho)+1) );
end

% Get PSD of P4 signal
[Pss_P4,w_P4] = periodogram(s_p4); %,ones(1024,1),1024);
f_plot_P4 = w_P4 / (2*pi*t_samp);

% Get autocorrelation of P4
xc_p4 = xcorr(s_p4);

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
%{
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
%}

figure(3)
plot(f_plot_P4,10*log10(Pss_P4));
title('P4 PSD');
xlabel('Frequency (hz)');
ylabel('Power (dB)');
%axis([0e6 400e6 -40 30]);
grid on;

figure(4)
plot(s_p4_lp);

figure(5)
plot(xc_lfm);
title('ACS, Chirp');

figure(6)
plot(xc_p4);
title('ACS, P4');