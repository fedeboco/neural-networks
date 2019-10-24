%Federico Pérez Boco 96777 - 2do cuat. 2018
%actualizar_asincronico(W, patron, iteraciones, plot_figure)
%calc_energia_umbral(s,W,umbral)
close all
clear all

N = 200; %Numero de neuronas
m = 10; %Numero de patrones
umbral = 0.7; %Umbral

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
        s(indices(i)) = signo(W(indices(i),:)*s - umbral);
        s_alterada(indices(i)) = signo(W_alterada(indices(i),:)*s_alterada - umbral);
        E = [E, calc_energia_umbral(s,W,umbral)];
        E_alterada = [E_alterada, calc_energia_umbral(s_alterada,W_alterada,umbral)];
    end   

end

plot(E)
hold on
plot(E_alterada)

legend('W completa','W alterada')
ylabel('Energía')
xlabel('Iteraciones')