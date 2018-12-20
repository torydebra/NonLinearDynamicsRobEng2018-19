function eqPoint = fnReazStability(a, b)
%% function fnReazStability(a, b)
% Function to analize the system fnReaz
%
% Input:
%   a = parameter alfa
%   b = parameter beta

% Find null isoclines
syms X Y Z

f1 = 5 * Z.^2 - 2 * X.^2 - 10 * X * Y;
f2 = a * Z - 0.1 * Y - 10 * X * Y;
f3 = - 0.0675 * Z - b * 0.0675 * (1 - X - Y - Z);

eqns(1) = f1 == 0;
eqns(2) = f2 == 0;
eqns(3) = f3 == 0;

S = solve(eqns, [X Y Z], 'Real', true);

nEqPoint = length(S.X);

eqPoint = zeros(nEqPoint,3);

for i=1:nEqPoint
    eqPoint(i,:) = real(double([S.X(i), S.Y(i), S.Z(i)]));
end

df1dX = diff(f1, X);
df1dY = diff(f1, Y);
df1dZ = diff(f1, Z);
df2dX = diff(f2, X);
df2dY = diff(f2, Y);
df2dZ = diff(f2, Z);
df3dX = diff(f3, X);
df3dY = diff(f3, Y);
df3dZ = diff(f3, Z);

J = [df1dX, df1dY, df1dZ; df2dX, df2dY, df2dZ; df3dX, df3dY, df3dZ];

Jpe = zeros(3, 3, nEqPoint);

eigValues = zeros(nEqPoint, 3);

for i=1:nEqPoint
    Jpe(:,:,i) = subs(J, [X Y Z], eqPoint(i,:));
    eigValues(i,:) = eig(Jpe(:,:,i));
    plotc(eigValues(i,:));
end
eqPoint
eigValues
