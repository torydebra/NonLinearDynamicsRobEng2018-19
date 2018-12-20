function [value, isterminal, direction] = maxReaz(t, x, a, b)
% Funzione 'event' per rilevare i massimi di x2 nel sistema
% a, b sono i parametri e sono *necessari* anche se non usati

% funzione di test: in questo caso corrisponde alla sezione di Poincare' 
% (punti in cui si annulla la derivata di x2)
value = - 0.0675 * x(3) - b * 0.0675 * (1 - x(1) - x(2) - x(3));

% l'integrazione non termina quando rileva l'evento
isterminal = 0;

% rileva l'evento solo in una direzione (prendo solo i massimi di x2 e non 
% i minimi)
direction = -1;