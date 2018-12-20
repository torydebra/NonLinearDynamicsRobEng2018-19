%VEDERE CHE SUCCE PER VARI PAR
close all
clearvars

params = [2.1, 2.1];

eqPoint = fnReazStability(params(1), params(2));

X0 = [0.2, 0.5 ,0.4];

opt = odeset('RelTol', 1.0e-5, ...
             'Refine', 32);
[T, X] = ode45(@fnReaz, [0 1000], X0, opt, params(1), params(2));
figure
plot3(X(1:1000,1),X(1:1000,2),X(1:1000,3),'g');
hold on
plot3(X(:,1),X(:,2),X(:,3),'b');
hold on
plot3(X(end,1),X(end,2),X(end,3),'ro');
xlabel('x');
ylabel('y');
zlabel('z');

% Grafico di x(t) nel tempo
figure
plot(T, X(:,1:3));
legend('x', 'y', 'z');