function h=plotFloq(x)
% x e' un vettore di numeri complessi
% disegna nel piano complesso dei pallini:
% verdi se la parte reale<0
% rossi se la parte reale>=0
re=real(x(:)); % re contiene tutte le parti reali di x
im=imag(x(:)); % im contiene tutte le parti immaginarie di x
hold on;
plot(re(re>=1),im(re>=1),'ro')
plot(re(re<1),im(re<1),'go')

% Traccia una linea nera in corrispondenza dell'asse immaginario
% Traccia una linea magenta in corrispondenza dell'asse reale
plot([0 0],[-2 2],'k',[-2 2],[0 0],'m');
xlabel('Re');
ylabel('Im');
axis tight;
hold off;
