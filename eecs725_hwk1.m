% EECS725 Homework 1

clear;
close all;
hold off;

% Universal constants
c = 3e8; % m/s

% Problem constraints
f = 1e9; % hz
l_x = 50 / 100; % m
l_y = 10 / 100; % m
v_x = 90; % m/s
h = 1000; % m
offset = -2000; % m
rcs = 1; % m^2

x = linspace(-2000,2000,100); % m

%% Problem 1
R = sqrt(h^2 + offset^2 + x.^2); % m

figure(1);
plot(x,R);
grid on;
axis([-2000 2000 0 4000]);
title('Radial distance vs. x');
xlabel('x (m)');
ylabel('R (m)');

%% Problem 2
phi = atan2(2000,-x); % radians

figure(2);
plot(x,phi*180/pi);
grid on;
axis([-2000 2000 0 180]);
title('Angle vs. x');
xlabel('x (m)');
ylabel('phi (degrees)');

%% Problem 3
theta = (45*pi/180) + atan2(-h,-offset); % radians
theta = theta * ones(size(x)); % vectorize

phi = phi - pi/2; % rotate phi 90 degrees to align to antenna vector
lambda = c/f; % m
beta_xz = lambda / l_x; % radians
beta_yz = lambda / l_y; % radians
beta_phi = beta_xz;
beta_theta = beta_yz;
Aeff = l_x * l_y; % m^2
Go = ((4*pi)/(lambda^2))*Aeff;
%G = Go * (sin(2.773.*theta/beta_theta)./(2.773.*theta/beta_theta)).^2 ...
%       .* (sin(2.773.*phi/beta_phi)./(2.773.*phi/beta_phi)).^2;
G = zeros(1,length(x));
for i=1:length(x)
    G(i) = Go * (sin(2.773*theta(i)/beta_theta)/(2.773*theta(i)/beta_theta))^2 ...
              * (sin(2.773*phi(i)/beta_phi)/(2.773*phi(i)/beta_phi))^2;
end

figure(3);
plot(x,10*log10(G));
grid on;
title('Antenna gain vs. x');
xlabel('x (m)');
ylabel('gain (dB)');

%% Problem 4
PrPt = (G.^2)*(lambda^2)*rcs ./ ( (4*pi)^3 * R.^4 );

figure(4);
plot(x,10*log10(PrPt));
grid on;
title('Power vs. x');
xlabel('x (m)');
ylabel('Power (dB)');

%% Problem 5
% Compute radial velocity
v_r = v_x * x./R; % m/s

% Compute Doppler shift
f_d = f * -v_r/c;

figure(5)
plot(x,f_d);
grid on;
title('Doppler shift vs. x');
xlabel('x (m)');
ylabel('Doppler shift (hz)');


