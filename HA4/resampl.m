function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

N = size(W,2);

W = W/sum(W);   %normalize

cum = cumsum(W);
samples = rand(N, 1);

j = zeros(1,N);

for i = 1:N
    j(i) = find(samples(i) <= cum, 1, 'first');
end

Xr = X(:,j);
Wr = ones(1,N)/N;


end
