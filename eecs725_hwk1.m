% EECS725 Homework 1

clear;
close all;
hold off;

% Universal constants
c = 3e8;        % speed of light (m/s)

% Problem constraints
f = 1e9;        % carrier frequency (hz)
l_x = 50 / 100; % antenna array width (m)
l_y = 10 / 100; % antenna array height (m)
v_x = 90;       % x speed (m/s)
h = 1000;       % elevation (m)
offset = -2000; % target offset (m)
sigma = 1;      % target rcs (m^2)

x = linspace(-2000,2000,100);   % position vectors (m)

%% Problem 1
% Compute radial distance
R = sqrt(h^2 + offset^2 + x.^2); % radial distance (m)

figure(1);
plot(x,R);
grid on;
axis([-2000 2000 0 3000]);
title('Radial distance vs. x');
xlabel('x (m)');
ylabel('R (m)');

%% Problem 2
% Compute offset angle
phi = atan2(2000,-x); % offset angle (radians)

figure(2);
plot(x,phi*180/pi);
grid on;
axis([-2000 2000 0 180]);
title('Offset angle vs. x');
xlabel('x (m)');
ylabel('phi (degrees)');

%% Problem 3
% Compute target elevation (theta) relative to dish orientation
theta = acos(h./R) - (45*pi/180); % radians

phi = phi - pi/2;       % rotate phi -90 degrees to align to antenna vector
lambda = c/f;           % wavelength (m)
beta_xz = lambda / l_x; % 3dB BW (radians)
beta_yz = lambda / l_y; % 3dB BW (radians)
beta_phi = beta_xz;
beta_theta = beta_yz;
Aeff = l_x * l_y;       % effective antenna array area (m^2)
Go = ((4*pi)/(lambda^2))*Aeff; % max antenna gain
G = Go * (sin(2.773.*theta/beta_theta)./(2.773.*theta/beta_theta)).^2 ...
      .* (sin(2.773.*phi  /beta_phi)  ./(2.773.*phi  /beta_phi))  .^2;

figure(3);
plot(x,10*log10(G));
grid on;
title('Antenna gain vs. x');
xlabel('x (m)');
ylabel('gain (dB)');

%% Problem 4
% Compute power ratio of received to transmit power from target reflection
PrPt = (G.^2)*(lambda^2)*sigma ./ ( (4*pi)^3 * R.^4 );

figure(4);
plot(x,10*log10(PrPt));
grid on;
title('Power ratio (Pr/Pt) vs. x');
xlabel('x (m)');
ylabel('Pr/Pt (dB)');

%% Problem 5
% Compute radial velocity
v_r = v_x * x./R; % radial velocity (m/s)

% Compute Doppler shift
f_d = f * -v_r/c; % doppler shift (hz)

figure(5)
plot(x,f_d);
grid on;
title('Doppler shift vs. x');
xlabel('x (m)');
ylabel('Doppler shift (hz)');


