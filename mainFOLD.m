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

eqPoint = fnReazStability(params(1), params(2))



%% Ricerca FOLD (backward of p1)
disp('PASSO 1: Ricerca della biforcazione di FOLD del punto di equilibrio P0=[-0.1169  -0.0019  0.0005]');
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto [1.6, 0.4]');
disp('             e si diminuisce il parametro p1')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(3,:),params,1)');

[x0,v0]=init_EP_EP(@reaz,eqPoint(3,:)',params',1);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',600); 
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward', 0);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp(['COMMENTO: E'' stata trovata una biforcazione FOLD (LP) per p1 = ', num2str(xE(4,sE(2).index))] );
disp('          (p2 e'' fissato a 0.4)');

figure;
cpl(xE,vE,sE,[4 1]);
xlabel('p_{1}')
ylabel('x_{1}')
zlabel('x_{2}')
title('Continuazione dell''equilibrio P0=[-0.1169  -0.0019  0.0005]')
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% Continuazione della biforcazione fold

% indice, nel vettore delle soluzioni, della soluzione corrispondente alla
% fold
indice = sE(2).index;
% coordinate dell'equilibrio nel punto di biforcazione
ptLP = xE(1:3,indice);
% valore del parametro nel punto di biforcazione
params(1) = xE(4,indice);
% inizializzo la continuazione della curva fold
[x0, v0] = init_LP_LP(@reaz, ptLP, params', [1,2]);
% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt,'Backward', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 400);
[xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt);
disp(' ');
disp(['COMMENTO: E'' stata trovata una biforcazione BT (BT) per p1,p2 = ', num2str(xLP(4,sLP(2).index)), ...
                       ' ', num2str(xLP(5,sLP(2).index))] );

% PLOTS per la CUSP

figure; hold on; view(3); grid on; box on;
sz = size(xLP(1,:));
minX = min(xLP(1,:))-0.5;
maxX = max(xLP(1,:))+0.5;
minP1 = min(xLP(4,:))-0.5;
maxP1 = max(xLP(4,:))+0.25;
minP2 = min(xLP(5,:))-0.25;
maxP2 = max(xLP(5,:))+0.25;

 
cpl(xLP,vLP,sLP,[4,5,1]);
plot3(xLP(4,:),xLP(5,:),repmat(minX,sz),'r','LineWidth',2); %plot projection on p1p2 plane


plot3(repmat(maxP1,sz),xLP(5,:),xLP(1,:),'r','LineWidth',2); %plot projection on p1X plane


plot3(xLP(4,:),repmat(maxP2,sz),xLP(1,:),'r','LineWidth',2); %plot projection on p2X plane

title('Continuazione della biforcazione FOLD al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2'); zlabel('x');
axis([minP1,maxP1,minP2,maxP2,minX,maxX]);

