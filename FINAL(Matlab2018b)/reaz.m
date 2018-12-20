function out = reaz
%
% Odefile of for the reaction
%
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
%
% -------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0; 10];
y0 = zeros(2,1);
handles = feval(@reaz);
options = odeset('Jacobian', handles(3), 'JacobianP',handles(4));
%
% -------------------------------------------------------------------------
           
function dfdt = fun_eval(t, x, a, b)
%
f1 = 5*x(3).^2-2*x(1).^2-10*x(1)*x(2);
f2 = a*x(3)-0.1*x(2)-10*x(1)*x(2);
f3 = -0.0675*x(3)-b*0.0675*(1-x(1)-x(2)-x(3));
dfdt = [f1; f2; f3];
%
% -------------------------------------------------------------------------

function jac = jacobian(t, x, a, b)
%
df1dXYZ = [-4*x(1)-10*x(2), -10*x(1), 10*x(3)];
df2dXYZ = [-10*x(2), -0.1-10*x(1), a];
df3dXYZ = [b*0.0675, b*0.0675, b*0.0675-0.0675];
jac = [df1dXYZ; df2dXYZ; df3dXYZ];
%
% -------------------------------------------------------------------------

function jacp = jacobianp(t, x, a, b)
%
jacp = [0, 0; x(3), 0; 0, -0.0675*(1-x(1)-x(2)-x(3))];
