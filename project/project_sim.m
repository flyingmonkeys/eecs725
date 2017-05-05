% EECS725 Project - Scenario simulation
% This simulation will compare performance of a linear frequency modulated
% (LFM), or "chirp", radar pulse transmission with a polyphase modulated P4
% transmission.

clear;
close all;

% Universal constants
c         = 3e8;       % speed of light (m/s)

% Problem constraints
f_start   = 10e9-(75e6/2); % chirp start frequency (hz)
f_end     = 10e9+(75e6/2); % chirp end frequency (hz)
Nc        = 128;           % pulse compression ratio = BW * tau

% Simulation parameters
samp_per_chip = 4;      % (samples/chip)
f_max         = 25;     % Range of Doppler frequencies to plot (normalized to f*tau)

A_t1      = 0.1;        % amplitude of return, target 1
A_t2      = 1.0;        % amplitude of return, target 2
A_t3      = 0.25;       % amplitude of return, target 3
d_t1      = 30;         % target 1 distance (m)
d_t2      = 90;         % target 2 distance (m)
d_t3      = 94;         % target 3 distance (m)
s_t1      = 4412;       % speed, target 1 (m/s)
s_t2      = 0;          % speed, target 2 (m/s)
s_t3      = 4412;       % speed, target 3 (m/s)

% Derived parameters
bw        = f_end - f_start;        % bandwidth (hz)
t_b       = 1 / bw;                 % P4 chip period (s)
fc        = (f_start + f_end) / 2;  % center (carrier) frequency for P4 (hz)
tau       = Nc / bw;                % pulse width (s)
t_samp    = t_b / samp_per_chip;    % simulation sampling period (s)
fd_max    = f_max / tau;            % Range of Doppler frequencies to plot (hz)
fd_t1     = 2 * fc * s_t1 / c;      % Doppler shift (hz), target 1
fd_t2     = 2 * fc * s_t2 / c;      % Doppler shift (hz), target 2
fd_t3     = 2 * fc * s_t3 / c;      % Doppler shift (hz), target 3

i_t1      = floor(2*d_t1/(c*t_samp));     % index of first target
i_t2      = floor(2*d_t2/(c*t_samp));     % index of second target
i_t3      = floor(2*d_t3/(c*t_samp));     % index of third target
t_burst   = linspace(0,tau,tau/t_samp)'; % vector of timestamps (s)

% Printout system parameters
fprintf(1,'Bandwidth = %5.2d Mhz\n',bw/1e6);
fprintf(1,'Pulse duration = %5.2d us\n',tau*1e6);
fprintf(1,'Range resolution = %5.2d m\n',c/(2*bw));
fprintf(1,'Simulation sampling frequency = %5.2d Mhz\n',1/(1e6*t_samp));


%% Chirp signal generation and analysis--------------------------------

% Generate chirp burst
k         = (f_end-f_start) / tau;           % chirp rate (hz/s)
theta_lfm = 2*pi*0.5*k.*(t_burst.^2);        % chirp phases
s_lfm     = cos( 2*pi*(f_start.*t_burst + 0.5*k.*(t_burst.^2)) );
s_lfm_bb  = exp( j*theta_lfm );              % baseband chirp pulse

% Get passband PSD of chirp
[Pss_lfm,w] = periodogram(s_lfm,[],1024);
f_plot  = w / (2*pi*t_samp);

% Get autocorrelation of chirp
xc_lfm_bb = abs(xcorr(s_lfm_bb,'biased'));

% Construct ambiguity plot (LFM baseband)
amb_lfm_bb      = zeros((Nc*2)-1,100);
doppler_lag_idx = 1;
for f_shift=0:fd_max/99:fd_max
    s_lfm_bb_shifted = exp( -j*2*pi*f_shift.*t_burst ) .* s_lfm_bb;
    xc = abs(xcorr(s_lfm_bb,s_lfm_bb_shifted,'biased'));
    amb_lfm_bb(:,doppler_lag_idx) = xc(length(s_lfm_bb)-Nc+1:length(s_lfm_bb)+Nc-1);
    doppler_lag_idx = doppler_lag_idx + 1;
end

%% P4 signal generation and analysis------------------------------------

% Generate P4 chip sequence
theta_p4 = zeros(Nc,1);
for i=2:Nc
    theta_p4(i) = pi * ( (((i-1)^2)/Nc) - (i-1) ); % P4 phases
end

% Generate P4 signal (baseband and passband signal)
s_p4_bb = zeros(length(t_burst),1);
for i=1:length(t_burst)
    s_p4_bb(i) = exp( j * theta_p4(mod(floor(i/samp_per_chip),Nc)+1) );
end
s_p4 = real( s_p4_bb .* exp(j*2*pi*fc.*t_burst) );

% Get passband PSD of P4 signal
[Pss_P4,w_P4] = periodogram(s_p4,[],1024);
f_plot_P4     = w_P4 / (2*pi*t_samp);

% Get autocorrelation of P4 (baseband)
xc_p4_bb = abs(xcorr(s_p4_bb,'biased'));

% Construct ambiguity plot (P4 baseband)
amb_p4_bb       = zeros((Nc*2)-1,100);
s_p4_bb_shifted = zeros(length(t_burst),1);
doppler_lag_idx = 1;
for f_shift=0:fd_max/99:fd_max
    s_p4_bb_shifted = exp( -j*2*pi*f_shift.*t_burst ) .* s_p4_bb;
    xc = abs(xcorr(s_p4_bb,s_p4_bb_shifted,'biased'));
    amb_p4_bb(:,doppler_lag_idx) = xc(length(s_p4_bb)-Nc+1:length(s_p4_bb)+Nc-1);
    doppler_lag_idx = doppler_lag_idx + 1;
end

%% Simulation start------------------------------------------------------

% Generate composite radar return signal (LFM), scenario 1
s_lfm_bb_tx = s_lfm_bb;
r = zeros(2*length(s_lfm_bb_tx),1);
r_t1 = r;
r_t2 = r;
r_t3 = r;
r_t1(i_t1:length(s_lfm_bb_tx)+i_t1-1) = A_t1 * s_lfm_bb_tx;
r_t2(i_t2:length(s_lfm_bb_tx)+i_t2-1) = A_t2 * s_lfm_bb_tx;
r_t3(i_t3:length(s_lfm_bb_tx)+i_t3-1) = A_t3 * s_lfm_bb_tx;
r = r_t1 + r_t2 + r_t3;

% Cross-correlate transmit and receive signals
xc_out_lfm  = flipud(abs(xcorr(s_lfm_bb_tx,r(1:length(s_lfm_bb_tx)),'biased')));
mixer_out_lfm_1 = xc_out_lfm(floor(length(xc_out_lfm)/2):length(xc_out_lfm));

% Generate composite radar return signal (LFM), scenario 2
r = zeros(2*length(s_lfm_bb_tx),1);
r_t1 = r;
r_t2 = r;
r_t3 = r;
r_t1(i_t1:length(s_lfm_bb_tx)+i_t1-1) = A_t1 * s_lfm_bb_tx .* exp( -j*2*pi*fd_t1.*t_burst );
r_t2(i_t2:length(s_lfm_bb_tx)+i_t2-1) = A_t2 * s_lfm_bb_tx .* exp( -j*2*pi*fd_t2.*t_burst );
r_t3(i_t3:length(s_lfm_bb_tx)+i_t3-1) = A_t3 * s_lfm_bb_tx .* exp( -j*2*pi*fd_t3.*t_burst );
r = r_t1 + r_t2 + r_t3;

% Cross-correlate transmit and receive signals
xc_out_lfm  = flipud(abs(xcorr(s_lfm_bb_tx,r(1:length(s_lfm_bb_tx)),'biased')));
mixer_out_lfm_2 = xc_out_lfm(floor(length(xc_out_lfm)/2):length(xc_out_lfm));

% Generate composite radar return signal (P4), scenario 1
s_p4_bb_tx = s_p4_bb;
r = zeros(2*length(s_p4_bb_tx),1);
r_t1 = r;
r_t2 = r;
r_t3 = r;
r_t1(i_t1:length(s_p4_bb_tx)+i_t1-1) = A_t1 * s_p4_bb_tx;
r_t2(i_t2:length(s_p4_bb_tx)+i_t2-1) = A_t2 * s_p4_bb_tx;
r_t3(i_t3:length(s_p4_bb_tx)+i_t3-1) = A_t3 * s_p4_bb_tx;
r = r_t1 + r_t2 + r_t3;

% Cross-correlate transmit and receive signals
xc_out_p4  = flipud(abs(xcorr(s_p4_bb_tx,r(1:length(s_p4_bb_tx)),'biased')));
mixer_out_p4_1 = xc_out_p4(floor(length(xc_out_p4)/2):length(xc_out_p4));


% Generate composite radar return signal (P4), scenario 2
r = zeros(2*length(s_p4_bb_tx),1);
r_t1 = r;
r_t2 = r;
r_t3 = r;
r_t1(i_t1:length(s_p4_bb_tx)+i_t1-1) = A_t1 * s_p4_bb_tx .* exp( -j*2*pi*fd_t1.*t_burst );
r_t2(i_t2:length(s_p4_bb_tx)+i_t2-1) = A_t2 * s_p4_bb_tx .* exp( -j*2*pi*fd_t2.*t_burst );
r_t3(i_t3:length(s_p4_bb_tx)+i_t3-1) = A_t3 * s_p4_bb_tx .* exp( -j*2*pi*fd_t3.*t_burst );
r = r_t1 + r_t2 + r_t3;

% Cross-correlate transmit and receive signals
xc_out_p4  = flipud(abs(xcorr(s_p4_bb_tx,r(1:length(s_p4_bb_tx)),'biased')));
mixer_out_p4_2 = xc_out_p4(floor(length(xc_out_p4)/2):length(xc_out_p4));

%% Plots------------------------------------------------------------

figure(1)
plot(linspace(0,length(mixer_out_p4_1)*c*t_samp/2,length(mixer_out_p4_1)),abs(mixer_out_p4_1));
hold on;
plot(linspace(0,length(mixer_out_p4_2)*c*t_samp/2,length(mixer_out_p4_1)),abs(mixer_out_p4_2));
grid on;
title('Matched filter response, P4, scenario 1 and 2');
xlabel('Distance (m)');
legend('Stationary','With Doppler');
axis([0 120 0 0.8]);


figure(2)
plot(linspace(0,length(mixer_out_lfm_1)*c*t_samp/2,length(mixer_out_lfm_1)),abs(mixer_out_lfm_1));
hold on;
plot(linspace(0,length(mixer_out_lfm_2)*c*t_samp/2,length(mixer_out_lfm_2)),abs(mixer_out_lfm_2));
grid on;
title('Matched filter response, LFM, scenario 1 and 2');
xlabel('Distance (m)');
legend('Stationary','With Doppler');
axis([0 120 0 0.8]);

figure(3)
xa = linspace(0,f_max,100);
ya = linspace(-0.5/samp_per_chip,0.5/samp_per_chip,Nc*2-1);
surf(xa,ya,amb_p4_bb);
title('Ambuguity plot, P4 baseband');
ylabel('lag (t/\tau)');
xlabel('Doppler shift (f_d\tau)');
grid on;

figure(4)
xa = linspace(0,f_max,100);
ya = linspace(-0.5/samp_per_chip,0.5/samp_per_chip,Nc*2-1);
surf(xa,ya,amb_lfm_bb);
title('Ambuguity plot, LFM baseband');
ylabel('lag (t/\tau)');
xlabel('Doppler shift (f_d\tau)');
grid on;

figure(5)
plot(20*log10(xc_lfm_bb));
hold on
plot(20*log10(xc_p4_bb));
title('ACS, LFM vs. P4')
xlabel('Lag');
ylabel('Power (dB)');
legend('LFM','P4');
grid on;

figure(6)
plot(f_plot_P4,20*log10(Pss_P4));
hold on
plot(f_plot,20*log10(Pss_lfm));
title('PSD of LFM and P4, Fc=200Mhz, BW=75Mhz');
xlabel('Frequency (hz)');
ylabel('Power (dB)');
legend('P4','LFM');
grid on;
