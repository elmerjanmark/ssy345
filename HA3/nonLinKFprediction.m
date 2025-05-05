function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
    n = size(x,1);
    
    switch type
        case 'EKF'
            [fx, Fx] = f(x);
            x = fx;
            P = Fx*P*Fx' + Q;
            
        case 'UKF'
            [SP,W] = sigmaPoints(x, P, type);
            x_est = zeros(n,1);
            P_est = Q;

            for i = 1:(2*n+1)
                [fx,~] = f(SP(:,i));
                
                fx_(:,i) = fx;
                
                x_est = x_est + fx*W(i);
                
            end
            
            for i = 1:(2*n+1)
                P_est = P_est + (fx_(:,i) - x_est)*(fx_(:,i) - x_est)'*W(i);
            end
            
            x = x_est;
            P = P_est;
            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
                
        case 'CKF'
            [SP,W] = sigmaPoints(x, P, type);
            
            x_est = zeros(n,1);
            P_est = Q;
            
            for i = 1:2*n
                [fx, ~] = f(SP(:,i));
                x_est = x_est + fx*W(i);
                fx_(:,i) = fx;
            end
            
            for i = 1:2*n
                P_est = P_est + (fx_(:,i)-x_est)*(fx_(:,i)-x_est)'*W(i);
            end
            
            x = x_est;
            P = P_est;
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end