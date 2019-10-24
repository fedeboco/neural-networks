%Federico Pérez Boco 96777 - 2do cuat. 2018
%actualizar_asincronico(W, patron, iteraciones, plot_figure)
close all
clear all

N = 200; %Numero de neuronas
m = 10; %Numero de patrones

%% Patrones aleatorios
P = signo(randn(N,m)); %Memorias

%% Aplico algoritmo de Hopfield
W = (P * P') - m * eye(N); %Enseño patrones a la red
W_alterada = W.*(rand(N) > (0.50)); %con 50% de conexiones faltantes

%% Genero estado inicial aleatorio
si = signo(randn(N,1)); %estado iniciales

%% Calculo energía
nro_actualizaciones = 3; %veces que actualizo cada neurona

s = si;
s_alterada = si;
E(1)=-1/2*s'*W*s;
E_alterada(1)=-1/2*s'*W_alterada*s;

for n = 1:nro_actualizaciones
    indices = randperm(N);
    for i = 1:N
        s(indices(i)) = signo(W(indices(i),:)*s);
        s_alterada(indices(i)) = signo(W_alterada(indices(i),:)*s_alterada);
        E = [E, -1/2*s'*W*s];
        E_alterada = [E_alterada, -1/2*s_alterada'*W_alterada*s_alterada];
    end   

end

plot(E)
hold on
plot(E_alterada)

legend('W completa','W alterada')
ylabel('Energía')
xlabel('Iteraciones')