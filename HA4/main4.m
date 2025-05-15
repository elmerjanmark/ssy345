%% 

close all; clear all; clc

rng(1);

% Sampling period 
T = 0.1; 
% Length of time sequence 
K = 800; 
% Allocate memory 
omega = zeros(1,K+1); 
% Set turn−rate at turns 
omega(200:350) = pi/301/T; 
omega(450:600) = pi/301/T; 
% Initial state 
x0 = [0 0 20 -pi/2 0]'; 
P0 = diag([10^2, 10^2, 10^2, (5*pi/180)^2, (5*pi/180)^2]);
% Allocate memory 
X = zeros(length(x0),K+1); 
X(:,1) = x0;
% Create true track 

for i=2:K+1 
    % Simulate 
    X(:,i) = coordinatedTurnMotion(X(:,i-1), T); 
    % Set turn−rate 
    X(5,i) = omega(i); 
end

sigma_v = 0.001;
sigma_w = 0.3*pi/180;
Q = diag([0, 0, sigma_v^2, sigma_w^2]);

sigma_phi1 = 0.5*pi/180;
sigma_phi2 = 0.5*pi/180;
R = diag([(sigma_phi1)^2, (sigma_phi2)^2]);

s1 = [160, -300]';
s2 = [420, -300]';

f = @(x) coordinatedTurnMotion(x,T);
h = @(x) dualBearingMeasurement(x,s1,s2);
Y = genNonLinearMeasurementSequence(X,h,R);

x_mes = zeros(1,K);
y_mes = zeros(1,K);
for i = 1:K
    [x_mes(i), y_mes(i)] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
end

[x_UKF, P_UKF, ~, ~] = nonLinearKalmanFilter(Y, x0, P0, f, Q, h, R, 'UKF');


figure
hold on; grid on;
plot(X(1,:), X(2,:), 'LineWidth', 1.5)
scatter(x_mes, y_mes, 'X')
scatter([s1(1) s2(1)],[s1(2) s2(2)],'k', 'filled')
plot(x_UKF(1,:), x_UKF(2,:), 'LineWidth', 1.5)
xlim([-100, 700])
ylim([-700 100])