function [x, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);

x_prev = x_0;
P_prev = P_0;

for k = 1:N

    % Prediction
    [x_pred, P_pred] = linearPrediction(x_prev, P_prev, A, Q);

    % Update
    y_k = Y(:,k);
    if (isnan(y_k(1))) && length(H)>1
        H_k = H(2,:);
        R_k = R(2,2);
        y_k = y_k(2);
    else
        H_k = H;
        R_k = R;
    end
    [x_upd, P_upd] = linearUpdate(x_pred, P_pred, y_k, H_k, R_k);

    x(:,k) = x_upd;
    P(:,:,k) = P_upd;

    % Prepare for next iteration
    x_prev = x_upd;
    P_prev = P_upd;
end