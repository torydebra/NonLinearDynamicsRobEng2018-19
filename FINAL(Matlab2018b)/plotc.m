function plotc(x)
% x e' un vettore di numeri complessi
% disegna nel piano complesso dei pallini:
% verdi se la parte reale<0
% rossi se la parte reale>=0
re=real(x(:)); % re contiene tutte le parti reali di x
im=imag(x(:)); % im contiene tutte le parti immaginarie di x
figure;
hold on;
plot(re(re>=0),im(re>=0),'ro')
plot(re(re<0),im(re<0),'go')

% Traccia una linea rossa in corrispondenza dell'asse immaginario
if min(im) < max(im)
    y=[min(im) max(im)];
else
    y=[-1 1];
end
plot([0 0],y,'r','linewidth',2);
xlabel('Re');
ylabel('Im');
