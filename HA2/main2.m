%% Scenario 1

%% a
rng(6)
Q = 1.5;
R = 3;

N = 35;

A = 1;
H = 1;

x_0 = 0;
mu = 2;
P_0 = 8;

X = genLinearStateSequence(x_0, P_0, A, Q, N);
Y = genLinearMeasurementSequence(X, H,R);

figure
plot(X)
hold on
plot(1:N,Y)
legend('State sequence X', 'Measurements Y')
xlabel('N')


%% b

[x, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

sigma_x = sqrt(squeeze(P)');
upper_bound = x + 3*sigma_x;
lower_bound = x - 3*sigma_x;

figure
hold on
plot(X, 'b')                          % True state
plot(Y, ':r')                           % Measurements
plot(x, 'k')         % Kalman estimate
plot(upper_bound, '--g')              
plot(lower_bound, '--g')              
legend('True state', 'Measurements', 'Estimate', '+3\sigma', '-3\sigma')
title('Kalman Filter Estimate with ±3σ Bound')
xlabel('Time step k')
ylabel('State')
grid on


k_vals = [1, 2, 4,30];
x_vals = linspace(-10, 10, 1000);

figure
hold on
colors = ['r', 'g', 'b', 'm'];
for i = 1:length(k_vals)
    k = k_vals(i);
    err_mean = x(k) - X(k);          % Estimation error
    err_std = sqrt(P(:,:,k));        % Standard deviation at time k
    pdf_vals = normpdf(x_vals, 0, err_std);

    plot(x_vals, pdf_vals, 'Color', colors(i), 'DisplayName', ['k = ' num2str(k)])
end

title('Error Densities at Selected Time Steps')
xlabel('Error')
ylabel('Density')
legend
grid on



%% c

x0_wrong = 12;
[x_wrong, P_wrong] = kalmanFilter(Y, x0_wrong, P_0, A, Q, H, R);

figure
hold on
plot(X, '--k')                          % True state
plot(x, '-r', 'LineWidth', 1)         % Kalman estimate
plot(x_wrong, '-b', 'LineWidth', 1)         % Kalman estimate         
legend('True state', 'Correct prior Estimate', 'Incorrect prior est')
title('Kalman Filter Estimate with incorrect prior mean')
xlabel('Time step k')
ylabel('State')
grid on


%% d

k = 15;

x_values = linspace(-5, 15, 1000);

% prior
prio = normpdf(x_values, x(k-1), sigma_x(k-1));

%prediction
[x_p, P_p] = linearPrediction(x(k-1), sigma_x(k-1)^2, A, Q);
pred = normpdf(x_values, x_p, sqrt(P_p));

% measurement
meas = normpdf(x_values, Y(k), R);

%update
[x_u, P_u] = linearUpdate(x_p, P_p, Y(k), H, R);
update = normpdf(x_values, x_u, sqrt(P_u));

figure
plot(x_values, prio)
hold on
grid on
plot(x_values, pred)
plot(x_values, meas)
plot(x_values, update)
legend('prior', 'prediction', 'measurements', 'update')
title('PDFs for kalman')


%% e



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scenario 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = load('SensorMeasurements.mat');

v0 = data.CalibrationSequenceVelocity_v0;
v10 = data.CalibrationSequenceVelocity_v10;
v20 = data.CalibrationSequenceVelocity_v20;

% Known velocities
true_v0 = 0;
true_v10 = 10;
true_v20 = 20;

C_v10 = mean(v10) / true_v10;
C_v20 = mean(v20) / true_v20;

C = mean([C_v10, C_v20])


%residuals
rv_v0 = v0/C -  true_v0;
rv_v10 = v10/C - true_v10;
rv_v20 = v20/C -  true_v20;

var_rv_v0 = var(rv_v0);
var_rv_v10 = var(rv_v10);
var_rv_v20 = var(rv_v20);

var_rv = mean([var_rv_v0, var_rv_v10, var_rv_v20])


figure;

subplot(1, 2, 1);
plot(v0, 'DisplayName', 'Measured (v=0)');
hold on;
plot(C * true_v0 * ones(size(v0)), 'r--', 'DisplayName', 'True (scaled)');
title('Velocity Measurements (v=0)');
xlabel('Sample');
ylabel('Velocity (m/s)');
legend;
hold off;

subplot(1, 2, 2);
plot(v10, 'DisplayName', 'Measured (v=10)');
hold on;
plot(v20, 'DisplayName', 'Measured (v=20)');
plot(C * true_v10 * ones(size(v10)), 'r--', 'DisplayName', 'True (scaled)');
plot(C * true_v20 * ones(size(v20)), 'g--', 'DisplayName', 'True (scaled)');
title('Velocity Measurements (v=10, 20)');
xlabel('Sample');
ylabel('Velocity (m/s)');
legend;


C_est = C;  % estimated C

noise10 = mvnrnd(0, var_rv_v10, 2000)';
noise20 = mvnrnd(0, var_rv_v20, 2000)';

% generate using model yv_k = C(v_k + rv_k)
fakeTrainSpeed10 = C_est * (10 + noise10);
fakeTrainSpeed20 = C_est * (20 + noise20);

% Plot and compare
figure
subplot(2,1,1)
plot(fakeTrainSpeed20, 'DisplayName', 'Generated (v=20)')
hold on
plot(v20, 'DisplayName', 'Measured (v=20)')
legend
xlabel('Sample')
ylabel('Speed (m/s)')
grid on

subplot(2,1,2)
plot(fakeTrainSpeed10, 'DisplayName', 'Generated (v=10)')
hold on
plot(v10, 'DisplayName', 'Measured (v=10)')
legend
xlabel('Sample')
ylabel('Speed (m/s)')
grid on


%% 2b

Y_seq = Generate_y_seq;
Y_seq(2,:) = Y_seq(2,:)/C; %scale the velocity

%% Constant velocity model C

T = 0.1;
var_p = 1;
var_rv = var_rv;
q_cv = 0.05;

A_cv = [1 T; 0 1];

H = [1  0;
     0  1];

Q_cv = 0.01*[T^3/3 T^2/2; T^2/2 T];

R_cv = diag([var_p, var_rv]);


x_0 = [0; 0];
P_0_cv = diag([var_p, var_rv]);

[x_c, P_c] = kalmanFilter(Y_seq, x_0, P_0_cv, A_cv, Q_cv, H, R_cv);


figure
plot(x_c(1, :), LineWidth=3)
hold on
grid on
plot(find(~isnan(Y_seq(1, :))), Y_seq(1, ~isnan(Y_seq(1, :))), '--k', LineWidth=3);
legend('estimated position', 'measured position')

figure
plot(x_c(2, :))
hold on
grid on
plot(Y_seq(2,:))
legend('estimated velocity', 'measured scaled velocity')


%% CA
% q_ca = 0.5;
x_0_ca = [0;0;0];

A_ca = [1 T 0.5*T^2;
        0 1 T;
        0 0 1];

% Measurement Matrix (measures position and velocity, not acceleration)
H_ca = [1 0 0;
        0 1 0];

% Measurement Noise Covariance (same sensor noise as for CV)
R_ca = diag([var_p, var_rv]);

Q_ca = 0.01*[T^5/20  T^4/8  T^3/6;
            T^4/8   T^3/3  T^2/2;
            T^3/6   T^2/2  T];

P_0_ca = diag([var_p, var_rv, 1.0]);

% Run CA Kalman Filter
[x_ca, P_ca] = kalmanFilter(Y_seq, x_0_ca, P_0_ca, A_ca, Q_ca, H_ca, R_ca);

figure
plot(x_ca(1, :), LineWidth=3)
hold on
grid on
plot(find(~isnan(Y_seq(1, :))), Y_seq(1, ~isnan(Y_seq(1, :))), '--k', LineWidth=3);
legend('estimated position', 'measured position')

figure
plot(Y_seq(2,:),'r')
hold on
grid on
plot(x_ca(2, :),'b', LineWidth=2)
legend( 'measured scaled velocity','estimated velocity')