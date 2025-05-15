function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Your code here
px = x(1);
py = x(2);

dx1 = px - s1(1);
dy1 = py - s1(2);

dx2 = px - s2(1);
dy2 = py - s2(2);

hx = [atan2(dy1, dx1);
      atan2(dy2, dx2)];

n = length(x);     
Hx = zeros(2, n);

r1_sq = dx1^2 + dy1^2;
r2_sq = dx2^2 + dy2^2;

Hx(1, 1) = -dy1 / r1_sq;   % d(theta1)/d(px)
Hx(1, 2) =  dx1 / r1_sq;   % d(theta1)/d(py)

Hx(2, 1) = -dy2 / r2_sq;   % d(theta2)/d(px)
Hx(2, 2) =  dx2 / r2_sq;   % d(theta2)/d(py)
end