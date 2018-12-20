function xdot = fnReaz(t, x, a, b)
%% function xdot = fnReaz(t, x, a, b)
% Chimical reaction model
%
% Input:
%   t = time instant (not used here but compulsory for odesolver)
%   x = state of the system (x, y, z)
%   a = parameter alfa
%   b = parameter beta

X = x(1,:);
Y = x(2,:);
Z = x(3,:);

n = size(x,2);

% Odesolver wants column vector
xdot = zeros(3,n);
xdot(1,:) = 5 * Z.^2 - 2 * X.^2 - 10 * X * Y;
xdot(2,:) = a * Z - 0.1 * Y - 10 * X * Y;
xdot(3,:) = - 0.0675 * Z - b * 0.0675 * (1 - X - Y - Z);
