syms x y
f1 = x*(3-x-2*y)
f2 = y*(2-x-y)


c = linspace(-10,10);
figure
plot(subs(solve(f1, x),y,c),c)
hold on
plot(c,subs(solve(f2,y),x,c))


x = sym('x', [1,3]);
syms a b

f1 = 5*x(3).^2 - 2*x(1).^2 - 10*x(1)*x(2);
f2 = a*x(3) - 0.1*x(2) - 10*x(1)*x(2);
f3 = - 0.0675*x(3) - b*0.0675*(1-x(1)-x(2)-x(3));
dfdt = [f1; f2; f3];

dfdt = subs(dfdt, [a,b], [1.4707, 0.4]);

equations(1) =  dfdt(1) == 0;
equations(2) =  dfdt(2) == 0;
equations(3) =  dfdt(3) == 0;

S = solve(equations, [x(1) x(2) x(3)], 'Real', true);

eq(1,:) = double(S.x1);
eq(2,:) = double(S.x2);
eq(3,:) = double(S.x3);



%%non serve in sta funza
%tspan = [0, 100];
%opt = odeset;
%[t, x] = ode45(@fnReaz,tspan,x0,opt,a,b);


