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

params = [2.1, 1.4];

eqPoint = fnReazStability(params(1), params(2));

eqPointChoosen = 2; % punto di equilibrio stabile

%% ricerca la biforcazione (hopf)
disp(['PASSO 1: Ricerca della biforcazione hopf del punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),')']);
disp('             e si aumenta il parametro p2')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(2,:),params,2)');

[x0,v0]=init_EP_EP(@reaz,eqPoint(eqPointChoosen,:)',params',2);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',0);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

figure;
cpl(xE,vE,sE,[4 1 2]);
xlabel('p_{2}')
ylabel('x_{1}')
zlabel('x_{2}')
title(['Continuazione dell''equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))])
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% continua il ciclo originato dalla Hopf rispetto a p2
disp(['PASSO 2: Continuazione del ciclo generato dalla biforcazione di Hopf del punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),') localizzato al PASSO 1']);
disp('             e si aumenta il parametro p2')
disp('premi un tasto...'); pause;

x1=xE(1:3,sE(2).index);
params(2)=xE(4,sE(2).index);

disp(' ');
disp('Calcolo del ciclo iniziale --> [x0,v0]=init_H_LC@reaz,x1,params,1,0.0001,20,4);')

[x0,v0]=init_H_LC(@reaz,x1,params',2,0.0001,20,4);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'MaxNumPoints',150);
opt=contset(opt,'Backward',0);

disp('Continuazione del ciclo --> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt)')
disp(' ');

[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);

%   grafico
figure;
plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);
title('Continuazione del ciclo al variare di $p_2$.','Interpreter','latex');
xlabel('$p_2$','Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');
zlabel('$x_2$','Interpreter','latex');

% Periodo del ciclo al variare del parametro.
figure;
plot(xlc(end,:),xlc(end-1,:),'k','LineWidth',2);
grid on;
xlabel('p_2'); ylabel('T');
title('Periodo del ciclo');