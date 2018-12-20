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

params = [1.6, 1];

eqPoint = fnReazStability(params(1), params(2))
%close all


%% ricerca la Hopf (punto eq3 nato dalla fold in p1=1.3)
%disp(['PASSO 1: Ricerca della biforcazione di Hopf super-critica del punto di equilibrio P3=', ...
   % num2str(eqPoint(3,1)), ' ' ,  num2str(eqPoint(3,2)), ' ',  num2str(eqPoint(3,3))] );
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto [1.4, 1.65]');
disp('             e si aumenta il parametro p1')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(1,:),params,1)');
[x0,v0]=init_EP_EP(@reaz,eqPoint(3,:)',params',2);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',600); 
opt=contset(opt,'Singularities',4);
opt = contset(opt,'Backward', 0);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp('COMMENTO: E'' stata trovata una biforcazione di Hopf per p1 = 1.56');
disp('          (p2 e'' fissato a 1.65)');

figure;
cpl(xE,vE,sE,[4 1]);
xlabel('p_{2}')
ylabel('x_{1}')
zlabel('x_{2}')
title(['Continuazione dell''equilibrio P3=', ...
    num2str(eqPoint(3,1)), ' ' ,  num2str(eqPoint(3,2)), ' ',  num2str(eqPoint(3,3))])
grid;