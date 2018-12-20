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

params = [1.6, 0.4];

eqPoint = fnReazStability(params(1), params(2));

eqPointChoosen = 2;

%% ricerca la biforcazione (fold)
disp(['PASSO 1: Ricerca della biforcazione fold del punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),')']);
disp('             e si diminuisce il parametro p1')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(1,:),params,1)');

[x0,v0]=init_EP_EP(@reaz,eqPoint(eqPointChoosen,:)',params',1);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',1);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp(['COMMENTO: E'' stata trovata una biforcazione fold per p1 = ',num2str(xE(4,sE(2).index))]);
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

%% Continuazione della biforcazione fold
% indice, nel vettore delle soluzioni, della soluzioni corrispondente alla
% fold
indice = sE(2).index;
% coordinate dell'equilibrio nel punto di biforcazione
ptLP = xE(1:3,indice);
% valore del parametro nel punto di biforcazione
params(1) = xE(4,indice);

disp(['PASSO 2: Continuazione della curva di biforcazione fold relativa al punto di equilibrio P0 = ',num2str(eqPoint(eqPointChoosen,:))]);
disp(['DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (',num2str(params(1)),',',num2str(params(2)),')','localizzato al PASSO 1']);
disp('             e si lasciano liberi entrambi i parametri p1 e p2')
disp('premi un tasto...'); pause;
disp(' ');

disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0] = init_LP_LP(@reaz, ptLP, params, [1,2])');
% inizializzo la continuazione della curva fold
[x0, v0] = init_LP_LP(@reaz, ptLP, params', [1,2]);

disp('Continuazione della biforcazione fold --> [xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt)')
disp(' ');
% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 200);
%opt = contset(opt,'Backward',1);
[xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt);

disp(' ');
disp(['COMMENTO: E'' stata trovata una biforcazione BT per (p1,p2) = (',num2str(xLP(4,sLP(2).index)),',',num2str(xLP(5,sLP(2).index)),')']);

figure; 
cpl(xLP,vLP,sLP,[4 5]);
title('Continuazione della biforcazione fold al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');