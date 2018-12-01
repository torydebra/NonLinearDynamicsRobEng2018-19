%% Poincare section
% Param to have a cycle
a = 1.9;
b = 1.6;


% integrazione numerica
x0 = [0.6147 0.148 0.575];
opt = odeset('RelTol', 1.0e-5, ...
             'Events', @maxReaz, ...
             'Refine', 32);
[T, X, Te, Xe, Ie] = ode45(@fnReaz, [0, 1000], x0, opt, a, b);

% piano (ydot = 0) che è la sezine di poincare
%[x,y]=meshgrid(linspace(-0.1,0.1,100),linspace(-0.1,0.1,100));
%z = (0.1/a) * y + (10/a)*x*y;

% piano (zdot = 0) che è la sezine di poincare
[x,y]=meshgrid(linspace(0,0.2,100),linspace(0.2,1.4,100));
z = (-b*(1-x-y))/(1-b);


% grafici
% Andamento nel tempo
figure;
plot(T,X(:,:));
legend('x', 'y', 'z');
hold on;
plot(Te,Xe(:,3),'go');
% Traiettoria completa nello spazio di stato, con intersezioni con la
% sezione di Poincare'
figure;
plot3(X(:,1),X(:,2),X(:,3),'b');
hold on;
h1=surface(x,y,z,'facecolor',[1.0000 0.5020 0],'edgecolor','none','facealpha',0.75);
plot3(Xe(:,1),Xe(:,2),Xe(:,3),'ro');
xlabel('x');
ylabel('y');
zlabel('z');


