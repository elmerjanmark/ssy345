function X = genNonLinearStateSequence(x_0, P_0, f, Q, N)
%GENNONLINEARSTATESEQUENCE generates an N+1-long sequence of states using a 
%    Gaussian prior and a nonlinear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

% Your code here

n = size(Q,1);
X = zeros(n, N+1);

% init state from gaussian prior
X(:, 1) = mvnrnd(x_0, P_0)';

for k = 2:N+1
    % eval motion model at current state
    [fx, ~] = f(X(:, k-1));
    
    %noise
    q = mvnrnd(zeros(n,1), Q)';
    
    % next state
    X(:, k) = fx + q;
end

end