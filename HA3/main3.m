%% 

close all; clear all; clc


% question 1
%% a)
close all; clear all; clc

rng(2)

x1 = [30; 30];
p1 = [(8/3)^2 0; 0 (4/3)^2];

x2 = [30; 5];
p2 = [(8/3)^2 0; 0 (1/3)^2];

% bearing sensors
s1 = [-25; 0];
s2 = [25; 0];

std = 0.1*pi/180;
R = [std^2 0; 0 std^2];

T = 1;
Ns_val = [100, 10000];

f = @(x) coordinatedTurnMotion(x,T);
h = @(x) dualBearingMeasurement(x,s1,s2);

for i = 1:2
    Ns = Ns_val(i);
    x1_s = mvnrnd(x1, p1, Ns)';
    x2_s = mvnrnd(x2, p2, Ns)';
    
    y1_s = genNonLinearMeasurementSequence(x1_s, h, R);
    y2_s = genNonLinearMeasurementSequence(x2_s, h, R);
    
    
    y1_mean = (1/Ns)*sum(y1_s')';
    y1_p    = (1/(Ns-1))*(y1_s - y1_mean)*(y1_s - y1_mean)';
    % y1_p = cov(y1_s')

    fprintf("y1_mean (Ns = %d): [%.6f; %.6f]\n", Ns, y1_mean(1), y1_mean(2));
    fprintf("y1_p (Ns = %d):\n%.6e  %.6e\n%.6e  %.6e\n\n", Ns, y1_p(1,1), y1_p(1,2), y1_p(2,1), y1_p(2,2));

    y2_mean = (1/Ns)*sum(y2_s')';
    y2_p    = (1/(Ns-1))*((y2_s - y2_mean)*(y2_s - y2_mean)');
    % y2_p = cov(y2_s')

    
    fprintf("y2_mean (Ns = %d): [%.6f; %.6f]\n", Ns, y2_mean(1), y2_mean(2));
    fprintf("y2_p (Ns = %d):\n%.6e  %.6e\n%.6e  %.6e\n\n", Ns, y2_p(1,1), y2_p(1,2), y2_p(2,1), y2_p(2,2));

    disp('===================================================================================')
        
    figure()
    subplot(1,2,1); hold on; grid on;

    title(sprintf('Dual-Bearing Uncertainty (Ns = %d)', Ns));
    xlabel('Bearing 1 (rad)');
    ylabel('Bearing 2 (rad)');

    % Plot samples
    scatter(y1_s(1,:), y1_s(2,:), 5, 'filled', 'DisplayName', sprintf('Samples, Ns = %d', Ns));

    % Plot mean
    plot(y1_mean(1), y1_mean(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Mean');
    % Plot 3σ ellipse
    xy1 = sigmaEllipse2D(y1_mean, y1_p, 3, 100);
    plot(xy1(1,:), xy1(2,:), 'r-', 'LineWidth', 2, 'DisplayName', '3σ ellipse');

    legend show

    subplot(1,2,2); hold on; grid on;

    title(sprintf('Dual-Bearing Uncertainty (Ns = %d)', Ns));
    xlabel('Bearing 1 (rad)');
    ylabel('Bearing 2 (rad)');

    % Plot samples
    scatter(y2_s(1,:), y2_s(2,:), 5, 'filled', 'DisplayName', sprintf('Samples, Ns = %d', Ns));
    % Plot mean
    plot(y2_mean(1), y2_mean(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Mean');
    % Plot 3σ ellipse
    xy2 = sigmaEllipse2D(y2_mean, y2_p, 3, 100);
    plot(xy2(1,:), xy2(2,:), 'r-', 'LineWidth', 2, 'DisplayName', '3σ ellipse');

    legend show
end


%% b

Q=R;

type_ekf = 'EKF';
type_ukf = 'UKF';

% For the first state density (x1, p1)
[y1_mean_ekf, y1_p_ekf] = nonLinKFprediction(x1, p1, h, Q, type_ekf);
[y1_mean_ukf, y1_p_ukf] = nonLinKFprediction(x1, p1, h, Q, type_ukf);


fprintf("y1_mean_EKF (Ns = %d): [%.6f; %.6f]\n", Ns, y1_mean_ekf(1), y1_mean_ekf(2));
fprintf("y1_mean_UKF (Ns = %d): [%.6f; %.6f]\n", Ns, y1_mean_ukf(1), y1_mean_ukf(2));

fprintf("y1_p_EKF (Ns = %d):\n%.6e  %.6e\n%.6e  %.6e\n\n", Ns, y1_p_ekf(1,1), y1_p_ekf(1,2), y1_p_ekf(2,1), y1_p_ekf(2,2));
fprintf("y1_p_UKF (Ns = %d):\n%.6e  %.6e\n%.6e  %.6e\n\n", Ns, y1_p_ukf(1,1), y1_p_ukf(1,2), y1_p_ukf(2,1), y1_p_ukf(2,2));

disp('====================================================================')

% For the second state density (x2, p2)
[y2_mean_ekf, y2_p_ekf] = nonLinKFprediction(x2, p2, h, Q, type_ekf);
[y2_mean_ukf, y2_p_ukf] = nonLinKFprediction(x2, p2, h, Q, type_ukf);

fprintf("y2_mean_EKF (Ns = %d): [%.6f; %.6f]\n", Ns, y2_mean_ekf(1), y2_mean_ekf(2));
fprintf("y2_mean_UKF (Ns = %d): [%.6f; %.6f]\n", Ns, y2_mean_ukf(1), y2_mean_ukf(2));

fprintf("y2_p_EKF (Ns = %d):\n%.6e  %.6e\n%.6e  %.6e\n\n", Ns, y2_p_ekf(1,1), y2_p_ekf(1,2), y2_p_ekf(2,1), y2_p_ekf(2,2));
fprintf("y2_p_UKF (Ns = %d):\n%.6e  %.6e\n%.6e  %.6e\n\n", Ns, y2_p_ukf(1,1), y2_p_ukf(1,2), y2_p_ukf(2,1), y2_p_ukf(2,2));



disp('====================================================================')



%% 
for i = 1:2
    if i == 1
        x_prior = x1;
        P_prior = p1;
        y_samples = y1_s;
        y_mean_mc = y1_mean; %monte carlo
        y_p_mc = y1_p;
        y_mean_ekf_val = y1_mean_ekf;
        y_p_ekf_val = y1_p_ekf;
        y_mean_ukf_val = y1_mean_ukf;
        y_p_ukf_val = y1_p_ukf;
        title_str = 'State Density 1';
    else
        x_prior = x2;
        P_prior = p2;
        y_samples = y2_s;
        y_mean_mc = y2_mean; %monte carlo
        y_p_mc = y2_p;
        y_mean_ekf_val = y2_mean_ekf;
        y_p_ekf_val = y2_p_ekf;
        y_mean_ukf_val = y2_mean_ukf;
        y_p_ukf_val = y2_p_ukf;
        title_str = 'State Density 2';
    end

    % EKF already done

    % UKF propagated sigma points
    n = length(x_prior);
    [SP_x_ukf, W_ukf] = sigmaPoints(x_prior, P_prior, 'UKF');
    SP_y_ukf = zeros(2, size(SP_x_ukf, 2)); % init

    for j = 1:size(SP_x_ukf, 2)
        SP_y_ukf(:, j) = h(SP_x_ukf(:, j));
    end

    figure();
    hold on; grid on;
    title(sprintf('Measurement - %s', title_str));
    xlabel('Bearing 1 (rad)');
    ylabel('Bearing 2 (rad)');

    % Plot samples of y
    scatter(y_samples(1,:), y_samples(2,:), 5, 'filled', 'DisplayName', 'Samples (y)');

    % Plot sample mean
    plot(y_mean_mc(1), y_mean_mc(2), 'kx', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'Mean (MC)');
    xy = sigmaEllipse2D(y_mean_mc, y_p_mc, 3, 100);
    plot(xy(1,:), xy(2,:), 'LineWidth', 2, 'DisplayName', 'SigmaEllipse samples');

    % Plot EKF approximated mean
    plot(y_mean_ekf_val(1), y_mean_ekf_val(2), 'mo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean (EKF)');
    xy = sigmaEllipse2D(y_mean_ekf_val, y_p_ekf_val, 3, 100);
    plot(xy(1,:), xy(2,:),':', 'LineWidth', 2, 'DisplayName', 'SigmaEllipse EKF');

    % Plot UKF approximated mean
    plot(y_mean_ukf_val(1), y_mean_ukf_val(2), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Mean (UKF)');
    xy = sigmaEllipse2D(y_mean_ukf_val, y_p_ukf_val, 3, 100);
    plot(xy(1,:), xy(2,:),'--','LineWidth', 2, 'DisplayName', 'SigmaEllipse UKF');

    % Plot propagated sigma points (UKF)
    scatter(SP_y_ukf(1,:), SP_y_ukf(2,:), 40, 'MarkerFaceColor','r', 'DisplayName', 'Sigma Points (UKF)');
    legend;
    hold off;
end



%% 2
clear all;
rng(2)
T = 1;
s1 = [-50; 0];
s2 = [50; 0];

f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);


sigma_v2 = 0.3^2;                  
sigma_omega2 = (0.3*pi/180)^2;    
sigma_phi2 = (0.1*pi/180)^2;       


Q = diag([0, 0, sigma_v2, 0, sigma_omega2]);

R = sigma_phi2 * eye(2);

% x0 = [0; 200; 2; 0; 0];
x0 = [200; 50; 2; 0; 0];
P0 = diag([(5/3)^2, (5/3)^2, 2^2, (0.1*pi/180)^2, (0.1*pi/180)^2]);

N = 100;

x_true = genNonLinearStateSequence(x0, P0, f, Q, N);
y_meas = genNonLinearMeasurementSequence(x_true, h, R);

[xf_EKF, Pf_EKF, ~, ~] = nonLinearKalmanFilter(y_meas, x0, P0, f, Q, h, R, 'EKF');
[xf_UKF, Pf_UKF, ~, ~] = nonLinearKalmanFilter(y_meas, x0, P0, f, Q, h, R, 'UKF');


for i = 1:length(y_meas)
    [x_mes(:,i), y_mes(:,i)] = getPosFromMeasurement(y_meas(1,i), y_meas(2,i), s1, s2);
end

figure; hold on; grid on;
title('True Path,Measurement, EKF & UKF Estimates with 3σ Ellipses');
xlabel('x [m]'); ylabel('y [m]');

% Plot true trajectory
plot(x_true(1,:), x_true(2,:), 'k', 'LineWidth',1.2);
plot(x_mes, y_mes, 'LineWidth', 1.2);


plot(xf_EKF(1,:), xf_EKF(2,:), '--', 'LineWidth', 1.2);

plot(xf_UKF(1,:), xf_UKF(2,:), '-.',  'LineWidth', 1.2);


scatter(s1(1), s1(2), 'filled')
scatter(s2(1), s2(2), 'filled')

for k = 10:10:N
    xy_ekf = sigmaEllipse2D(xf_EKF(1:2,k), Pf_EKF(1:2,1:2,k), 3, 50);
    plot(xy_ekf(1,:), xy_ekf(2,:), 'r', 'LineWidth', 1.5);

    xy_ukf = sigmaEllipse2D(xf_UKF(1:2,k), Pf_UKF(1:2,1:2,k), 3, 50);
    plot(xy_ukf(1,:), xy_ukf(2,:), 'b--', 'LineWidth', 1.5);
end

% Plot sensors


legend('True states', 'Measurement', 'EKF Estimate', 'UKF Estimate', 'sensor 1', 'sensor 2' , '3σ EKF', '3σ UKF')


%% c

x0_all = [[0; 200; 2; 0; 0],[200; 50; 2; 0; 0]];

for i = 1:2
    x0 = x0_all(:,i);
    err_x_EKF = []; 
    err_y_EKF = [];
    err_x_UKF = []; 
    err_y_UKF = [];

    for k = 1:100
        %generate
        X = genNonLinearStateSequence(x0, P0, f, Q, N);
        Y = genNonLinearMeasurementSequence(X, h, R);

        [xf_EKF,Pf_EKF,~,~] = nonLinearKalmanFilter(Y, x0, P0, f, Q, h, R, 'EKF');
        [xf_UKF,Pf_UKF,~,~] = nonLinearKalmanFilter(Y, x0, P0, f, Q, h, R, 'UKF');

        err_x_EKF = [err_x_EKF, xf_EKF(1,:) - X(1,2:end)];
        err_y_EKF = [err_y_EKF, xf_EKF(2,:) - X(2,2:end)];

        err_x_UKF = [err_x_UKF, xf_UKF(1,:) - X(1,2:end)];
        err_y_UKF = [err_y_UKF, xf_UKF(2,:) - X(2,2:end)];
    end

    mean_errx_EKF = mean(err_x_EKF);
    mean_erry_EKF = mean(err_y_EKF);
    mean_errx_UKF = mean(err_x_EKF);
    mean_erry_UKF = mean(err_y_EKF);

    sigma_errx_EKF = sqrt(var(err_x_EKF));
    sigma_erry_EKF = sqrt(var(err_y_EKF));
    sigma_errx_UKF = sqrt(var(err_x_UKF));
    sigma_erry_UKF = sqrt(var(err_y_UKF));
    
        
    figure; 
    titles = {'EKF x-error', 'EKF y-error', 'UKF x-error', 'UKF y-error'};
    errors = {err_x_EKF, err_y_EKF, err_x_UKF, err_y_UKF};

    for j = 1:4
        subplot(2,2,j); hold on; grid on;
        e = errors{j};
        histogram(e, 'Normalization', 'pdf');

        mu = mean(e);
        sigma = std(e);

        x_vals = linspace(mu - 4*sigma, mu + 4*sigma, 100);

        plot(x_vals, normpdf(x_vals, mu, sigma), 'r', 'LineWidth', 2);
        xline(mu, '--y', 'LineWidth', 1.5);
        title(sprintf('Case %d %s', i, titles{j}));
        legend('Calculated', 'Gaussian fit', 'Mean');
        xlim([mu - 4*sigma, mu + 4*sigma])
    end


end


%%  3


T = 0.1;
K = 800;
omega = zeros(1,K+1);
omega(200:350) = pi/301/T;
omega(450:600) = pi/301/T; 
% Initial state 
x0 = [0 0 20 -pi/2 0]'; 
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


