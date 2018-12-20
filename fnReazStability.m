function eqPoint = fnReazStability(a, b)
%% function fnReazStability(a, b) 
% function to analize the stability of system fnReaz
%
% Inputs:
%   a, b - parameters

if nargin < 2
     error("insert param");
end

%% Find null isoclines 
syms X Y Z;
f1 = 5*Z.^2 - 2*X.^2 - 10*X*Y;
f2 = a*Z - 0.1*Y - 10*X*Y;
f3 = - 0.0675*Z - b*0.0675*(1-X-Y-Z);
equations(1) =  f1 == 0;
equations(2) =  f2 == 0;
equations(3) =  f3 == 0;

S = solve(equations, [X Y Z], 'Real', true);  %only real solutions are eq points
nEqPoints = length(S.X);
eqPoint = zeros (nEqPoints, 3); 

for i = 1:nEqPoints
    eqPoint(i,:) = real(double([S.X(i), S.Y(i), S.Z(i)]));
end

df1dx = diff(f1, X);
df1dy = diff(f1, Y);
df1dz = diff(f1, Z);
df2dx = diff(f2, X);
df2dy = diff(f2, Y);
df2dz = diff(f2, Z);
df3dx = diff(f3, X);
df3dy = diff(f3, Y);
df3dz = diff(f3, Z);

J = [df1dx, df1dy, df1dz; df2dx, df2dy, df2dz; df3dx, df3dy, df3dz];

Jpe = zeros (3, 3, nEqPoints);
eigValues = zeros(nEqPoints, 3);
for i = 1:nEqPoints
    Jpe(:,:,i) = subs (J, [X, Y, Z], eqPoint(i,:));
    eigValues(i,:) = eig(Jpe(:,:,i));
    plotc(eigValues(i,:))
end
eigValues






