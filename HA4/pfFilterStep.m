function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
%PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
%
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k

N = size(X_kmin1, 2);       
n = size(X_kmin1, 1);      

X_pred = proc_f(X_kmin1);                               
process_noise = mvnrnd(zeros(n,1), proc_Q, N)';         
X_k = X_pred + process_noise;


%predicted measurements
y_hat = meas_h(X_k);

% likelyhood p(y_k | x_k^i)
W = mvnpdf(yk', y_hat', meas_R)';

W_k = W_kmin1 .* W;
W_k = W_k / sum(W_k);

end