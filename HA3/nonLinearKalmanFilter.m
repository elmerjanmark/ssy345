function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
%NONLINEARKALMANFILTER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%                       [fx,Fx]=f(x) 
%                       Takes as input x (state) 
%                       Returns fx and Fx, motion model and Jacobian evaluated at x
%   Q           [n x n] Process noise covariance
%   h                   Measurement model function handle
%                       [hx,Hx]=h(x,T) 
%                       Takes as input x (state), 
%                       Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%

% Your code here. If you have good code for the Kalman filter, you should re-use it here as
% much as possible.

N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
xf = zeros(n,N);
Pf = zeros(n,n,N);

xp = zeros(n,N);
Pp = zeros(n,n,N);

x_prev = x_0;
P_prev = P_0;

for k = 1:N

    % Prediction
    [x_pred, P_pred] = nonLinKFprediction(x_prev, P_prev, f, Q, type);
    xp(:,k) = x_pred;
    Pp(:,:,k) = P_pred;
    
    y_k = Y(:,k);
    [x_upd, P_upd] = nonLinKFupdate(x_pred, P_pred, y_k, h, R, type);

    xf(:,k) = x_upd;
    Pf(:,:,k) = P_upd;

    % Prepare for next iteration
    x_prev = x_upd;
    P_prev = P_upd;

end