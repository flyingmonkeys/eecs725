% EECS725 Project - Ambiguity Plots
% This simulation will compare ambiguity functions of a linear frequency
% modulated (LFM), or "chirp", radar pulse transmission with a polyphase
% modulated P4 pulse for a pulse compression ratio of 64.

clear;
close all;

% Problem constraints
f_start   = 100e6;      % chirp start frequency (hz)
f_end     = 175e6;      % chirp end frequency (hz)
Nc        = 64;         % pulse compression ratio = BW * tau

% Simulation parameters
samp_per_chip = 3;      % (samples/chip)
f_max         = 25;     % Range of Doppler frequencies to plot (normalized to f*tau)

% Derived parameters
bw        = f_end - f_start;        % bandwidth (hz)
t_b       = 1 / bw;                 % P4 chip period (s)
tau       = Nc / bw;                % pulse width (s)
t_samp    = t_b / samp_per_chip;    % simulation sampling period (s)
fd_max    = f_max / tau;            % Range of Doppler frequencies to plot (hz)

t_burst   = linspace(0,tau,tau/t_samp)'; % vector of timestamps (s)

%% Chirp signal generation and analysis--------------------------------

% Generate chirp burst
k         = (f_end-f_start) / tau;           % chirp rate (hz/s)
theta_lfm = 2*pi*0.5*k.*(t_burst.^2);        % chirp phases
s_lfm_bb  = exp( j*theta_lfm );              % baseband chirp pulse

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

% Generate P4 signal (baseband)
s_p4_bb = zeros(length(t_burst),1);
for i=1:length(t_burst)
    s_p4_bb(i) = exp( j * theta_p4(mod(floor(i/samp_per_chip),Nc)+1) );
end

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

%% Plots------------------------------------------------------------

figure(1)
xa = linspace(0,f_max,100);
ya = linspace(-0.5/samp_per_chip,0.5/samp_per_chip,Nc*2-1);
surf(xa,ya,amb_p4_bb);
title('Ambuguity plot, P4 baseband');
ylabel('lag (t/\tau)');
xlabel('Doppler shift (f_d\tau)');
zlabel('\mid\chi(f_d,\tau)\mid');
grid on;

figure(2)
xa = linspace(0,f_max,100);
ya = linspace(-0.5/samp_per_chip,0.5/samp_per_chip,Nc*2-1);
surf(xa,ya,amb_lfm_bb);
title('Ambuguity plot, LFM baseband');
ylabel('lag (t/\tau)');
xlabel('Doppler shift (f_d\tau)');
zlabel('\mid\chi(f_d,\tau)\mid');
grid on;

figure(3)
stairs(theta_p4);
title('P4 sequence');
xlabel('Chip number');
ylabel('Phase (radians)');
grid on;

figure(4)
plot(20*log10(xc_lfm_bb));
hold on
plot(20*log10(xc_p4_bb));
title('ACS, LFM vs. P4')
xlabel('Lag');
ylabel('Power (dB)');
legend('LFM','P4');
grid on;
