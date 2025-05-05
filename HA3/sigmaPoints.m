function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%


n = size(x,1);

    switch type        
        case 'UKF'
            SP = zeros(n, 2 * n + 1);
            SP(:, 1) = x;
            
            W = zeros(1, 2 * n + 1);
            W(1) = 1 - n/3;
            W(2:end) = (1 - W(1)) / (2 * n);
            
            Psqrt = chol(P, 'lower');
            
            for i = 1:n
                SP(:,i+1)     = x + sqrt(n/(1-W(1))) * Psqrt(:,i);
                SP(:,i+1+n)   = x - sqrt(n/(1-W(1))) * Psqrt(:,i);
            end
        case 'CKF'
            
            SP = zeros(n, 2 * n);
            W = zeros(1, 2 * n);

            W(:) = 1 / (2 * n);
            
            Psqrt = chol(P, 'lower');
            
            for i = 1:n
                SP(:,i)     = x + sqrt(n) * Psqrt(:,i);
                SP(:,i+n)   = x - sqrt(n) * Psqrt(:,i);
            end
            
        otherwise
            error('Incorrect type of sigma point')
    end

end