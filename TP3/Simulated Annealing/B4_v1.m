%Federico Pérez Boco 96777 - Redes Neuronales

clear all
close all
%Págs 116 y 120

%% Entrada

%XOR

%salida deseada
yd = [0 ; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0];

%patrón de entrada
p = [0 0 0 0;
     0 0 0 1;
     0 0 1 0;
     0 0 1 1;
     0 1 0 0;
     0 1 0 1;
     0 1 1 0;
     0 1 1 1;
     1 0 0 0;
     1 0 0 1;
     1 0 1 0;
     1 0 1 1;
     1 1 0 0;
     1 1 0 1;
     1 1 1 0;
     1 1 1 1];

 
%% Configuración

%Cte prob
kb = 1.38e-23;
%Número de entradas
Li = length(p(1,:));
%Número de patrones
Lp = length(p(:,1));
%Número de salidas
Lo = length(yd(1,:));
%Elementos de capas i-ésima (1 capa de entradas, k capas de neuronas, 1 capa de salidas)
N = [Li 3 Lo] +1;
Nt = length(N); %Número de elementos/Cs
%Beta función sigmoide
B = 1/2; %tipico 1/2 o 1
%Constante de aprendizaje nu
n = 0.1; %entre 0 y 1
%Inicializo pesos con valores pequeños
for i = 1:Nt-1
    Wo{i} = rand(N(i+1),N(i))./100;  
end
%Umbral de ECM
U = 0.01;
%Máxima cant de iteraciones
iter_MAX = 10000;
%Factor de amortiguamiento de cambio en matriz de pesos
fA = 10;
%Iteraciones por temperatura
iT = 100;

%% Procesamiento

iter = 1;
Ti = 5;
T = Ti;
deltaT = 0.01;
ECM_aux = inf;

while(T > 0.01 && iter < iter_MAX)
    
    indice = randperm(length(p(:,1))); %indice aleatorio para patrones
    
    %Itero varias veces para cada temperatura
    for cambios = 1:iT
        for i = 1:Nt-1
            Wn{i} = Wo{i} + randn(N(i+1),N(i))/fA;  
        end  
        
        ECM = 0;
        
        % Proceso todos los patrones %
        for j = 1:Lp

            % Proceso salida %
            V{1} = [p(indice(j),:)];
            for C = 1:Nt-1 %avanza en capas C=CAPA
                for k = 1:N(C+1)-1 %calcula para cada neurona excepto umbral
                    H{C}(k) = Wn{C}(k,:)*[V{C} 1]'; %suma para neurona k-esima, 1 por umbral
                end
                V{C+1} = sig2(H{C},B); %sigmoide en todas las neuronas de la capa
            end

            % Actualizo ECM para patrón j-ésimo y lo acumulo %
            ECM = ECM + sum((yd(indice(j),:)-[V{end}])*(yd(indice(j),:)-[V{end}])');

        end 

        %actualizo W----------------
        P = exp( (ECM_aux-ECM)/(kb*T) ); %empieza a variar mas con baja temp
        if(ECM <= ECM_aux)
            Wo = Wn;
            ECM_aux = ECM;
        elseif(P > rand)
            Wo = Wn;
            ECM_aux = ECM;
        end

    end
    ECM_aux
    T_grafico(iter) = T;
    ECM_grafico(iter) = ECM_aux;
    iter = iter + 1;
    T = T - deltaT*T/Ti %caida de la temperatura amortiguada

end

%% Ploteo ECM

figure
plot(T_grafico(2:end), log(ECM_grafico(2:end)))
title('ECM')
ylabel('log(ECM)')
xlabel('Temperatura')
grid on

%% Pruebo el algoritmo

for i=1:length(p)
    g(i) = calc_salida(Nt, p(i,1:end), Wn, N, B);
end

g