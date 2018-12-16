function [M, m] = dsfindpeaks(s)
% trova gli indici M dei massimi e m dei minimi di s

% Funziona cercando gli indici dove la derivata di s cambia segno

% Calcolo numericamente la derivata
ds = diff(s);

% Utilizzo solo il segno: 
% ds := 1 se ds>0
% ds := -1 se ds<0
ds = sign(ds);

% ds e' fatto in questo modo:
% 1 1 1 1 1 1             1 1 1 1 1 1
%
%             -1 -1 -1 -1
% allora rifaccio la differenza dell'i-esimo con il precedente
ds = diff(ds);

% ottengo un nuovo ds fatto cosi':
%                    2
% 0 0 0 0 0    0 0 0   0 0 0 0 0
%           -2
% a questo punto mi basta guardare il segno di ds per vedere in quali
% posizioni si trovano i minimi e i massimi
m = find(ds > 0) + 1;
M = find(ds < 0) + 1;