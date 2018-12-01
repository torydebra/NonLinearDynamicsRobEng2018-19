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

params = [1.4, 1.9]; % p2 1.65 original

eqPoint = fnReazStability(params(1), params(2))
close all


%% ricerca la Hopf (punto eq3 nato dalla fold in p1=1.3)
disp(['PASSO 1: Ricerca della biforcazione di Hopf super-critica del punto di equilibrio P3=', ...
    num2str(eqPoint(3,1)), ' ' ,  num2str(eqPoint(3,2)), ' ',  num2str(eqPoint(3,3))] );
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto [1.4, 1.65]');
disp('             e si aumenta il parametro p1')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(1,:),params,1)');
[x0,v0]=init_EP_EP(@reaz,eqPoint(3,:)',params',1);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',60); 
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward', 0);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp('COMMENTO: E'' stata trovata una biforcazione di Hopf per p1 = 1.56');
disp('          (p2 e'' fissato a 1.65)');

figure(3);
cpl(xE,vE,sE,[4 1]);
xlabel('p_{1}')
ylabel('x_{1}')
zlabel('x_{2}')
title(['Continuazione dell''equilibrio P3=', ...
    num2str(eqPoint(3,1)), ' ' ,  num2str(eqPoint(3,2)), ' ',  num2str(eqPoint(3,3))])
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% Continuazione della biforcazione Hopf 
%Trova una BT che è relativa alla seconda HOPF e LP trovato al passo
%precedente
disp(['Ricerca della curva di biforcazione di Hopf relativa al punto di equilibrio P3=', ...
    num2str(eqPoint(3,1)), ' ' ,  num2str(eqPoint(3,2)), ' ',  num2str(eqPoint(3,3))]);
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (1.56, 1.65)');
disp('             e si lasciano liberi entrambi i parametri p1 e p2')
disp('premi un tasto...'); pause;
disp(' ');

x1=xE(1:3,sE(2).index);
params(1)=xE(4,sE(2).index);

[x0,v0]=init_H_H(@reaz,x1,params',[1 2]);

[xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt);

disp(' ');
disp('COMMENTO: E stata trovata una BT per (p1,p2) = (2.22, 3.38)');

figure(4);
cpl(xH,vH,sH,[4 5]);
xlabel('p_{1}')
ylabel('p_{2}')
title('Continuazione della biforcazione di Hopf localizzata al PASSO 1')
grid;


disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% Continuazione Ciclo generato dalla HOPF in 1.56; 1.65




