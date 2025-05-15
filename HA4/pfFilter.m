function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Non-resampled Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.

    K = size(Y, 2);
    n = length(x_0);
    
    xfp = zeros(n,K);
    Pfp = zeros(n,n,K);
    
    Xp = zeros(n, N, K);
    Wp = zeros(N, K);
    
    X_k = mvnrnd(x_0, P_0, N)';
    W_k  = ones(1,N)/N;
    
    
    for k = 1:K
        
        [X_k, W_k] = pfFilterStep(X_k, W_k, Y(:,k), proc_f, proc_Q, meas_h, meas_R);
        
        if bResample
            [X_k, W_k, j] = resampl(X_k, W_k);
        end
        
        xfp(:, k) = X_k * W_k'; %mean
        Pfp(:, :, k) = (X_k - xfp(:, k)) * diag(W_k) * (X_k - xfp(:, k))'; %cov
    
        Xp(:, :, k) = X_k;
        Wp(:, k) = W_k';
        
        if ~isempty(plotFunc)
            plotFunc(X_k, W_k, k);
        end
        
    end
end