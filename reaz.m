%% file for matcont functions
function out = reaz
%
% Standard ode file of reaction 
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

function [tspan,y0,options] = init
tspan = [0; 10];
y0 = [0;0];
handles = feval(@reaz);
options = odeset('Jacobian',handles(3),'JacobianP',handles(4));

function dfdt = fun_eval(t,x,a,b)
f1 = 5*x(3).^2 - 2*x(1).^2 - 10*x(1)*x(2);
f2 = a*x(3) - 0.1*x(2) - 10*x(1)*x(2);
f3 = - 0.0675*x(3) - b*0.0675*(1-x(1)-x(2)-x(3));
dfdt = [f1; f2; f3];

function jac = jacobian(t,x,a,b)
df1dxyz = [-4*x(1)-10*x(2),    -10*x(1),    10*x(3)];
df2dxyz = [-10*x(2),    -0.1-10*x(1),    a];
df3dxyz = [b*0.0675,    b*0.0675,    -0.0675+b*0.0675];
jac = [df1dxyz; df2dxyz; df3dxyz];


function jacp = jacobianp(t,x,a,b)
df1dab = [0,    0];
df2dab = [x(3),    0];
df3dab = [0,    -0.0675*(1-x(1)-x(2)-x(3))];
jacp = [df1dab; df2dab; df3dab];

% function hess = hessians(t,x,a,b)
% hessf1 = [-4, -10, 0; -10, 0, 0; 0, 0, 10];
% hessf2 = [0, -10, 0; -10, 0, 0; 0, 0, 0];
% hessf3 = zeros(3);
% hess = cat(3,hessf1,hessf2,hessf3);

