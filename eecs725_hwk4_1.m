% EECS725 Homework 4 Problem 1

clear;
close all;

% Universal constants
c = 3e8;       % speed of light (m/s)

% Problem constraints
lambda = 2e-2;                          % radar wavelength (m)
tx_pos = [-23578e3 -20050e3 22254e3];   % transmitter position (m)
tx_vel = [0 0 0];                       % transmitter velocity (m/s)
rx_pos = [0 0 4.35e3];                  % receiver position (m)
rx_vel = [0 300 0];                     % receiver velocity (m/s)

% UNCOMMENT FOR EXAMPLE TEST CONSTRAINTS
%{
lambda = c/200e6;
tx_pos = [-6e3 -6e3 500];
rx_pos = [0 0 2e3];
rx_vel = [0 100 0]; % receiver velocity (m/s)
%}

% Simulation parameters
pos_samp = 10;       % simulation sampling spacing (m)

%% Start simulation

% Form planar surface
X = -6e3:pos_samp:6e3; % (m)
Y = -6e3:pos_samp:6e3; % (m)
[mx my] = meshgrid(X,Y); % (m)

% Compute range to planar surface
r0 = zeros(length(X),length(Y)); % range to target, time 0
for m=1:length(X)
    for n=1:length(Y)
        surface_pos = [X(m) Y(n) 0];
        r0(m,n) = norm(tx_pos - surface_pos) + norm(surface_pos - rx_pos);
    end
end
%
% Compute new range to planar surface based on aircraft movement
r1     = zeros(length(X),length(Y)); % range to target, time 1
dt     = pos_samp / norm(rx_vel); % elapsed time over the sampling space
rx_pos = rx_pos + (pos_samp*(rx_vel / norm(rx_vel))); % move the aircraft
for m=1:length(X)
    for n=1:length(Y)
        surface_pos = [X(m) Y(n) 0];
        r1(m,n) = norm(tx_pos - surface_pos) + norm(surface_pos - rx_pos);
    end
end

% Compute Isodoppler values
fb = -(1/lambda)*(r1-r0)/dt;

% Compute Isopower values
P = zeros(length(X),length(Y)); % Power values (dB)
for m=1:length(X)
    for n=1:length(Y)
        surface_pos = [X(m) Y(n) 0];
        P(m,n) = 20*log10( 1.0 / (norm(tx_pos - surface_pos) * norm(surface_pos - rx_pos)) );
    end
end

%% Plots-----------
figure(1)
contour(mx,my,r0'/1000,'LevelStep',0.74);
title('Isorange');
xlabel('X (m)');
ylabel('Y (m)');
grid on;
cb0 = colorbar;
cb0.Label.String = 'Range (km)';

figure(2)
contour(mx,my,fb'/1000,'LevelStep',1.5);
%contour(mx,my,fb','LevelStep',10);
title('Isodoppler');
xlabel('X (m)');
ylabel('Y (m)');
grid on;
cb1 = colorbar;
cb1.Label.String = 'Doppler Shift (khz)';

figure(3)
contour(mx,my,P','LevelStep',1);
title('Isopower');
xlabel('X (m)');
ylabel('Y (m)');
grid on;
cb2 = colorbar;
cb2.Label.String = 'Power variation (dB)';