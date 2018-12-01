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

params = [1.4707, 0.4];

eqPoint = fnReazStability(params(1), params(2));
close all


%% ricerca la Hopf
disp('PASSO 1: Ricerca della biforcazione di Hopf super-critica del punto di equilibrio P0=[0.64,-0.06,-0.28]');
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto [1.4707, 0.4]');
disp('             e si aumenta il parametro p1')
disp('premi un tasto...'); pause; disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(1,:),params,1)');

[x0,v0]=init_EP_EP(@reaz,eqPoint(1,:)',params',1);

disp('Continuazione dell''equilibrio --> [xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt)')
disp(' ');
opt=contset;

opt=contset(opt,'MaxNumPoints',20000);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward', 1);

[xE,vE,sE,hE,fE]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp('COMMENTO: E'' stata trovata una biforcazione di Hopf super-critica per p1 = 1400');
disp('          (p2 e'' fissato a 0.4)');

figure(3);
cpl(xE,vE,sE,[4 1 2]);
xlabel('p_{1}')
ylabel('x_{1}')
zlabel('x_{2}')
title('Continuazione dell''equilibrio P0=[0.64,-0.06,-0.28]')
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% continua la Hopf
disp('PASSO 2: Ricerca della curva di biforcazione di Hopf relativa al punto di equilibrio P0=[0.64,-0.06,-0.28]');
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (1400,0.4) localizzato al PASSO 1');
disp('             e si lasciano liberi entrambi i parametri p1 e p2')
disp('premi un tasto...'); pause;
disp(' ');

x1=xE(1:3,sE(2).index);
params(1)=xE(4,sE(2).index);

disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_H_H(@reaz,x1,params,[1 2])');

[x0,v0]=init_H_H(@reaz,x1,params',[1 2]);

disp('Continuazione della biforcazione di Hopf --> [xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt)')
disp(' ');

[xH,vH,sH,hH,fH]=cont(@hopf,x0,v0,opt);

disp(' ');
disp('COMMENTO: Sono state trovate due biforcazioni LP per (p1,p2) = (11.4,50.2) e per (p1,p2) = (-8.9,306.16)');

figure(4);
cpl(xH,vH,sH,[4 5]);
xlabel('p_{1}')
ylabel('p_{2}')
title('Continuazione della biforcazione di Hopf localizzata al PASSO 1')
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

%% continua il ciclo originato dalla Hopf rispetto a p1
disp('PASSO 3: Continuazione del ciclo stabile generato dalla biforcazione di Hopf del punto di equilibrio P0=[0.64,-0.06,-0.28]');
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (1400,0.4) localizzato al PASSO 1');
disp('             e si aumenta il parametro p1')
disp('premi un tasto...'); pause;

x1=xE(1:3,sE(2).index);
params(1)=xE(4,sE(2).index);

disp(' ');
disp('Calcolo del ciclo iniziale --> [x0,v0]=init_H_LC@reaz,x1,params,1,0.0001,20,4);')

[x0,v0]=init_H_LC(@reaz,x1,params',1,0.0001,20,4);
opt=contset(opt,'Singularities',0);
opt = contset(opt,'MaxNumPoints',100);

disp('Continuazione del ciclo --> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt)')
disp(' ');

[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);

disp(' ');
disp('COMMENTO: NON e'' stata trovata alcuna biforcazione del ciclo e l''analisi si interrompe per (p1,p2) = (?,0.4)');
disp('premi un tasto...'); pause;
disp(' ');

% figure(3);
% title('Continuazione del ciclo a partire dall''equilibrio P0=[0.64,-0.06,-0.28]')
% dim=[3,1]; %dimensione dello spazio di stato e numero di parametri attivi
% xdim=dim(1)*1;
% pdim=dim(2)+1;
% figure(3); hold on
% for i=1:size(xlc,2)
% plot3(repmat(xlc(end,i),size(xlc(1:xdim:end-pdim,i))), ...
%     xlc(1:xdim:end-pdim,i),xlc(2:xdim:end-pdim,i), 'b')
% end
% view(3);
% plot3(repmat(xlc(end,slc(2).index),size(xlc(1:xdim:end-pdim,slc(2).index))), ...
%     xlc(1:xdim:end-pdim,slc(2).index),xlc(2:xdim:end-pdim,slc(2).index), 'or');
% grid on;
% xlabel('p_{1}')
% ylabel('x_{1}')
% zlabel('x_{2}')

%   grafico
figure(3);
figure(3); hold on
plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);
hold on;
cpl(xE,vE,sE,[4 1 2]);
title('Continuazione del ciclo al variare di $p$.','Interpreter','latex');
xlabel('$p$','Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');
zlabel('$x_2$','Interpreter','latex');

disp('premi un tasto...'); pause;
disp(' ');

disp('PASSO 4: Ricerca della biforcazione di Hopf sub-critica del punto di equilibrio P0=[0.64,-0.06,-0.28]');
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (?,?)');
disp('             e si diminuisce il parametro p1')
disp('premi un tasto...'); pause;

params=[50, 60];

disp(' ');
disp('Calcolo del punto iniziale --> [x0,v0]=init_EP_EP(@reaz,eqPoint(1,:),params,1)')

[x0,v0]=init_EP_EP(@reaz,eqPoint(1,:)',params',1);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Backward',1);
opt=contset(opt,'MaxNumPoints',20000);

disp('Continuazione dell''equilibrio --> [xE1,vE1,sE1,hE1,fE1]=cont(@equilibrium,x0,[],opt)')
disp(' ');

[xE1,vE1,sE1,hE1,fE1]=cont(@equilibrium,x0,v0,opt);

disp(' ');
disp('COMMENTO: E'' stata trovata una biforcazione di Hopf sub-critica per p1 = 0');
disp('          (p2 e'' fissato a 1)');
disp('premi un tasto...'); pause;
disp(' ');

figure(5);
cpl(xE1,vE1,sE1,[4 1 2]);
xlabel('p_{1}')
ylabel('x_{1}')
zlabel('x_{2}')
title('Continuazione dell''equilibrio P0=[0.64,-0.06,-0.28]')
grid;

disp('premi un tasto...'); pause;
disp(' ');
disp(' ');

disp('PASSO 5: Continuazione del ciclo instabile generato dalla biforcazione di Hopf del punto di equilibrio P0=(0,0)');
disp('DESCRIZIONE: Nel piano dei parametri (p1,p2) si parte dal punto (0,1) localizzato al PASSO 4');
disp('             e si diminuisce il parametro p1')
disp('premi un tasto...'); pause;

x1=xE1(1:3,sE1(3).index);
params(1)=xE1(4,sE1(3).index);

disp(' ');
disp('Calcolo del ciclo iniziale --> [x0,v0]=init_H_LC(@bautin,x1,p,[1],0.0001,20,4)')

[x0,v0]=init_H_LC(@reaz,x1,params',1,0.0001,20,4);

opt=contset(opt,'MinStepsize',1.0e-9);
opt = contset(opt,'VarTolerance',1e-6);
opt = contset(opt,'FunTolerance',1e-6);
opt = contset(opt,'TestTolerance',1e-6);

opt = contset(opt,'MaxNumPoints',130);

disp('Continuazione del ciclo --> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt)')
disp(' ');

[xlc3,vlc3,slc3,hlc3,flc3]=cont(@limitcycle,x0,v0,opt);

disp(' ');
disp('COMMENTO: E'' stata trovata una biforcazione tangente di cicli per p1 = -0.25');
disp('          (p2 e'' fissato a 1)');
disp('premi un tasto...'); pause;
disp(' ');
