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

%% Cuspide
params = [1.4707, 0.4];

eqPoint = fnReazStability(params(1), params(2));

eqPointChoosen = 1;

[x0,v0]=init_EP_EP(@reaz,eqPoint(eqPointChoosen,:)',params',1);

opt=contset;
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',1);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

% indice, nel vettore delle soluzioni, della soluzioni corrispondente alla
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
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 150);
opt = contset(opt,'Backward',1);
[xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt);

figure(3); 
cpl(xLP,vLP,sLP,[4 5]);

%% Bogdanov-tackens 1 (fold)
params = [1.6, 0.4];

eqPoint = fnReazStability(params(1), params(2));

eqPointChoosen = 2;

[x0,v0]=init_EP_EP(@reaz,eqPoint(eqPointChoosen,:)',params',1);

opt=contset;

opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',1);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

% indice, nel vettore delle soluzioni, della soluzioni corrispondente alla
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
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 30);
[xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt);

figure(3); 
hold on;
cpl(xLP,vLP,sLP,[4 5]);
title('Continuazione della biforcazione fold al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');

% inizializzo la continuazione della curva fold
[x0, v0] = init_LP_LP(@reaz, ptLP, params', [1,2]);

% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 100);
opt = contset(opt, 'Backward', 1);
[xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt);

figure(3); 
hold on;
cpl(xLP,vLP,sLP,[4 5]);
title('Continuazione della biforcazione fold al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');


%% Bogdanov-tackens 1 (Hopf)
params = [1.4, 1.65];

eqPoint = fnReazStability(params(1), params(2));

eqPointChoosen = 2; % punto di equilibrio stabile

[x0,v0]=init_EP_EP(@reaz,eqPoint(eqPointChoosen,:)',params',1);

opt=contset;
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',0);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

% indice, nel vettore delle soluzioni, della soluzioni corrispondente alla
% fold
indice = sE(2).index;
% coordinate dell'equilibrio nel punto di biforcazione
ptH = xE(1:3,indice);
% valore del parametro nel punto di biforcazione
params(1) = xE(4,indice);

% inizializzo la continuazione della curva fold
[x0, v0] = init_H_H(@reaz, ptH, params', [1,2]);

% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 35);
opt = contset(opt,'Backward',0);
[xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt);

figure(3);
hold on;
cpl(xH,vH,sH,[4 5]);
title('Continuazione della biforcazione fold al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');

% inizializzo la continuazione della curva fold
[x0, v0] = init_H_H(@reaz, ptH, params', [1,2]);

% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 10);
opt = contset(opt,'Backward',1);
[xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt);

figure(3);
hold on;
cpl(xH,vH,sH,[4 5]);
title('Continuazione della biforcazione fold al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');

%% Bogdanov-tackens 2 (fold)
% indice, nel vettore delle soluzioni, della soluzioni corrispondente alla
% fold
indice = sE(4).index;
% coordinate dell'equilibrio nel punto di biforcazione
ptLP = xE(1:3,indice);
% valore del parametro nel punto di biforcazione
params(1) = xE(4,indice);

% inizializzo la continuazione della curva fold
[x0, v0] = init_LP_LP(@reaz, ptLP, params', [1,2]);

% opzioni di continuazione
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'MinStepSize', 1e-4);
opt = contset(opt, 'MaxNumPoints', 30);
opt = contset(opt,'Backward',1);
[xLP,vLP,sLP,hLP,fLP] = cont(@limitpoint, x0, v0, opt);

figure(3);
hold on;
cpl(xLP,vLP,sLP,[4 5]);
title('Continuazione delle biforcazioni al variare di (p_1,p_2).');
xlabel('p_1'); ylabel('p_2');
plot(-1.5, 0.1, 'g*');
plot(-1.4, 0.1, 'k*');
plot(-1.5, 0.5, 'g*');
plot(-1.4, 0.5, 'k*');
plot(-1.5, 0.6, 'g*');
plot(-1.4, 0.6, 'r*');
plot(0, 0.6, 'g*');
plot(0.1, 0.6, 'k*');
plot(2, 0.6, 'g*');
plot(2.1, 0.6, 'k*');
plot(2, 0.9, 'g*');
plot(2.1, 0.9, 'r*');
plot(1.8, 2.75, 'g*');
plot(1.9, 2.75, 'k*');
plot(1.8, 2.5, 'r*');
plot(1.9, 2.5, 'r*');
plot(2.8, 2.5, 'g*');
plot(2.9, 2.5, 'k*');
circle(1.564,1.65,0.07);
circle(1.594,1.65,0.08);
circle(1.641,1.65,0.09);
circle(1.675,1.65,0.1);