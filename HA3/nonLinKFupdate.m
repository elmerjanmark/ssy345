function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%
    n = size(x,1);
    switch type
        case 'EKF'
            [hx, Hx] = h(x);
            
            S = Hx*P*Hx' + R;
            K = P*Hx'*inv(S);
            
            P = P - K*S*K';
            x = x + K*(y - hx);
        case 'UKF'
            [SP,W] = sigmaPoints(x, P, type);
            
            y_est = 0;
            for i = 1:(2*n+1)
                [hx,~] = h(SP(:,i));
                hx_(:,i) = hx;
                
                y_est = y_est + hx*W(i);
            end
            
            P_xy = 0;
            S_est = R;
            
            for i = 1:(2*n+1)
                P_xy = P_xy + (SP(:,i) - x)*(hx_(:,i) - y_est)'*W(i);

                S_est = S_est + (hx_(:,i)-y_est) * (hx_(:,i)-y_est)' * W(i);
            end
            
            %update
            x = x + P_xy*inv(S_est)*(y - y_est);
            P = P - P_xy*inv(S_est)*P_xy';

            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
            
        case 'CKF'
            [SP,W] = sigmaPoints(x, P, type);
    
            y_est = 0;
            for i = 1:2*n
                [hx,~] = h(SP(:,i));
                hx_(:,i) = hx;
                
                y_est = y_est + hx*W(i);
            end
            
            P_xy = 0;
            S_est = R;
            
            for i = 1:(2*n)
                P_xy = P_xy + (SP(:,i) - x)*(hx_(:,i) - y_est)'*W(i);

                S_est = S_est + (hx_(:,i)-y_est) * (hx_(:,i)-y_est)' * W(i);
            end
            
            %update
            x = x + P_xy*inv(S_est)*(y - y_est);
            P = P - P_xy*inv(S_est)*P_xy';
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

