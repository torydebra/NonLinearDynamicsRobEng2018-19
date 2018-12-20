%%  Funzione per la lyapunov()
function f=fnReaz4Lyapuv(t,X,par)

a = par(1);
b = par(2);
x=X(1); y=X(2); z=X(3);

Y= [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];
f=zeros(9,1);
f(1)= 5*z.^2 - 2*x.^2 - 10*x*y;
f(2)= a*z - 0.1*y - 10*x*y;
f(3)= - 0.0675*z - b*0.0675*(1-x-y-z);

df1dxyz = [-4*x-10*y,    -10*x,    10*z];
df2dxyz = [-10*y,    -0.1-10*x,    a];
df3dxyz = [b*0.0675,    b*0.0675,    -0.0675+b*0.0675];
Jac = [df1dxyz; df2dxyz; df3dxyz];

f(4:12)=Jac*Y;