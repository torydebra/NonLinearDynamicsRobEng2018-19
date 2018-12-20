clearvars
close all
clc

addpath('Cl_ matcont4p2');
addpath('Cl_ matcont4p2/Equilibrium');
addpath('Cl_ matcont4p2/Continuer');
addpath('Cl_ matcont4p2/LimitPoint');
addpath('Cl_ matcont4p2/MultilinearForms');
addpath('Cl_ matcont4p2/Hopf');
addpath('Cl_ matcont4p2/LimitCycle');

params = [1.4, 1.65];

eqPoint = fnReazStability(params(1), params(2));

eqPointChoosen = 2; % punto di equilibrio stabile

%% ricerca la biforcazione (hopf)
disp(['PASSO 1: Ricerca della biforcazione hopf del punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),')']);
disp('             e si aumenta il parametro p1')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(2,:),params,1)');

[x0,v0]=init_EP_EP(@reaz,eqPoint(eqPointChoosen,:)',params',1);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',0);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp(['COMMENTO: Sono state trovate due biforcazioni hopf per p1 = ',num2str(xE(4,sE(2).index)),'e p1 = ',num2str(xE(4,sE(3).index))]);
disp(['          Inoltre e'' stata trovata una biforcazione fold per p1 = ',num2str(xE(4,sE(4).index))]);
disp(['          (p2 e'' fissato a ',num2str(params(2)),')']);

figure;
cpl(xE,vE,sE,[4 1 2]);
xlabel('p_{1}')
ylabel('x_{1}')
zlabel('x_{2}')
title(['Continuazione dell''equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))])
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% Continuazione di una delle due biforcazioni hopf
% indice, nel vettore delle soluzioni, della soluzioni corrispondente alla
% fold
indice = sE(2).index;
% coordinate dell'equilibrio nel punto di biforcazione
ptH = xE(1:3,indice);
% valore del parametro nel punto di biforcazione
params(1) = xE(4,indice);

disp(['PASSO 2: Continuazione della curva di biforcazione hopf relativa al punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),') localizzato al PASSO 1']);
disp('             e si lasciano liberi entrambi i parametri p1 e p2')
disp('premi un tasto...'); pause;
disp(' ');

disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0] = init_H_H(@reaz, ptLP, params, [1,2])');
% inizializzo la continuazione della curva fold
[x0, v0] = init_H_H(@reaz, ptH, params', [1,2]);

disp('Continuazione della biforcazione hopf --> [xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt)')
disp(' ');
% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 100);
opt = contset(opt,'Backward',1);
[xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt);

disp(' ');
disp(['COMMENTO: E'' stata trovata una biforcazione BT per (p1,p2) = (',num2str(xH(4,sH(2).index)),',',num2str(xH(5,sH(2).index)),')']);

figure;
cpl(xH,vH,sH,[4 5]);
title('Continuazione della biforcazione fold al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% continua il ciclo originato dalla Hopf rispetto a p1
disp(['PASSO 3: Continuazione del ciclo generato dalla biforcazione di Hopf del punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),') localizzato al PASSO 1']);
disp('             e si aumenta il parametro p1')
disp('premi un tasto...'); pause;

x1=xE(1:3,sE(2).index);
params(1)=xE(4,sE(2).index);

disp(' ');
disp('Calcolo del ciclo iniziale --> [x0,v0]=init_H_LC@reaz,x1,params,1,0.0001,20,4);')

[x0,v0]=init_H_LC(@reaz,x1,params',1,0.0001,20,4);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'MaxNumPoints',130);
opt=contset(opt,'Backward',0);

disp('Continuazione del ciclo --> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt)')
disp(' ');

[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);

disp(' ');
disp(['COMMENTO: NON e'' stata trovata alcuna biforcazione del ciclo e l''analisi si interrompe per (p1,p2) = (',num2str(params(1)),',',num2str(params(2)),')']);
disp('premi un tasto...'); pause;
disp(' ');

%   grafico
figure;
plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);
title('Continuazione del ciclo al variare di $p_1$.','Interpreter','latex');
xlabel('$p_1$','Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');
zlabel('$x_2$','Interpreter','latex');

% Rappresento i moltiplicatori di Floquet
figure;
title('moltiplicatori di Floquet del ciclo');
plotFloq(flc(end-2:end,10));
hold on;
% Cerchio di raggio unitario
[x y]=meshgrid(-2:0.1:2,-2:0.1:2);
contour(x,y,x.^2+y.^2-1,[0 0]);
axis([-2 2 -2 2]);
hold off;