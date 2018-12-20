%% Poincare section
% Param to have a cycle a 2.1
a = 2.1;
b = 1.7;

% integrazione numerica
x0 = [0.6147 0.148 0.575];
opt = odeset('RelTol', 1.0e-5, ...
             'Events', @maxReaz, ...
             'Refine', 32);
[T, X, Te, Xe, Ie] = ode45(@fnReaz, [0, 1000], x0, opt, a, b);

% piano (ydot = 0) che è la sezione di poincare
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

% Traiettoria a regime nello spazio di stato, con intersezioni con la
% sezione di Poincare'
figure;
tstart = 136000;
tend = floor(9*length(T)/10); % per levare as much as possible the transitorio
plot3(X(tend:end,1),X(tend:end,2),X(tend:end,3),'b');
hold on;
h1=surface(x,y,z,'facecolor',[1.0000 0.5020 0],'edgecolor','none','facealpha',0.75);
starte = floor(9*length(Te)/10);
plot3(Xe(starte:end,1),Xe(starte:end,2),Xe(starte:end,3),'ro');
xlabel('x');
ylabel('y');
zlabel('z');
figure
plot(Xe(starte:end,1),Xe(starte:end,3),'ro');
xlabel('x');
ylabel('z');

%% Diagramma biforcazione
% Crea diagramma di biforcazione al variare di r

% N e' il numero di valori del parametro di biforcazione che sara'
% utilizzato per costruire il diagramma.
N = 100;
% Raffino la tolleranza per l'integrazione. Calcolero' i massimi di x1
% usando non la funzionalita' Events delle ode, ma una funzione ad hoc,
% basata sul risultato dell'integrazione numerica.
opt = odeset('RelTol', 1.0e-5);

% Vettore dei valori di r per i quali calcolo il diagramma di biforcazione
A = linspace(1.91, 2.1, N);

figure;
axes('NextPlot', 'Add');
xlabel('a');
ylabel('z^*');
title('Diagramma di biforcazione');
for a=A
    disp(a)
    % Transitorio (da scartare)
    % Anche la scelta della durata del transitorio e' critica
    [T, X] = ode45(@fnReaz, [tstart tend], x0, opt, a, b);
    % Riassegno la condizione iniziale
    x0 = X(end, :);
    
    % Classifico solo questi dati
    [T, X] = ode45(@fnReaz, [0 20], x0, opt, a, b);
    
    % Ricavo numericamente i massimi di x1 a partire dai risultati
    % dell'integrazione numerica
    [p, t] = dsfindpeaks(X(:, 3));
    
    % Visualizzo i massimi di x1, ossia le intersezioni della traiettoria
    % con la sezione di Poincare' scelta, di equazione dx1/dt = 0
    % N.B.: ones e' *molto* piu' veloce di repmat
    % VERIFICA:
    %     tic; r(ones(size(p))); toc
    %     Elapsed time is 0.000082 seconds.
    %     tic; repmat(r, size(p)); toc
    %     Elapsed time is 0.026891 seconds.
    plot(a(ones(size(p))), X(p, 3), '.', 'MarkerSize', 20);    
end
% Disegno rette verticali in corrispondenza di 3 valori di r che
% corrispondono a comportamenti asintotici diversi
% plot([250 250], [40 50], 'r');
% plot([222 222], [40 50], 'g');
% plot([216.5 216.5], [40 50], 'c');

%% Lyapunov ESPO
a= 1.9;
b = 2;

% Dimensione del vettore di stato
n = 3;
% Passo sulla variabile t per la procedura di rinormalizzazione di
% Gram-Schmidt
stept = 0.5;
% Passo di stampa sulla finestra di comando
ioutp = 10;
% Chiamo la funzione che effettua l'integrazione numerica e il calcolo
% degli esponenti di Lyapunov dalla traiettoria calcolata
[T1,Res]=mylyapunov(n,@fnReaz4Lyapuv,@ode45,tstart,stept,tend,x0,ioutp,[a,b]);

% Grafico degli esponenti di Lyapunov al variare del tempo
figure;
plot(T1,Res);
hold on;
plot([tstart tend],[0,0],'k');
ylim([-10, 1]);



