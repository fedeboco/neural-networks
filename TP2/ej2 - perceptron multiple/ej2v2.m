%Federico Pérez Boco 96777 - Redes Neuronales

clear all
close all
%Págs 116 y 120

%% Entrada

%XOR

%salida deseada
yd = [0 1 1 0]';

%patrón de entrada
p = [0 0;
     0 1;
     1 0;
     1 1];

%umbral
p = [p ones(size(p(:,1)))];
 
%% Configuración
 
%Número de entradas
Li = length(p(1,:));

%Número de salidas
Lo = length(yd(1,:));

%Elementos de capas i-ésima (1 capa de entradas, k capas de neuronas, 1 capa de salida)
N = [Li 3 Lo];
Nt = length(N); %Número de elementos/Cs

%Beta función sigmoide
B = 1/2; %tipico 1/2 o 1

%Constante de aprendizaje nu
n = 0.7; %entre 0 y 1

%Inicializo pesos con valores pequeños
for i = 1:Nt-1
    W{i} = rand(N(i+1),N(i))./100;  
end

%Umbral de ECM
U = 0.005;

%% Procesamiento

O = zeros(1,4);
ECM(1) = inf;
iter = 1;

while(ECM(iter) > U)
    
    indice = randperm(length(p(:,1))); %indice aleatorio para patrones
    for j = 1:indice %permuta entre patrones

        %forward-----------------------------
        V{1} = p(indice(j),:);
        for C = 1:Nt-1 %avanza en capas C=CAPA
            
            for k = 1:N(C+1) %calcula para cada neurona
                H{C}(k) = W{C}(k,:)*V{C}';
            end
            V{C+1} = sig(H{C},B);
        end
        
        %backward----------------------------
        delta{Nt} = dsig(H{end},B)*(yd(indice(j))-V{end});
        for C = Nt-1:-1:2 %para cada capa excepto la última

             for k = 1:N(C) %capa actual
                 for l = 1:N(C+1) %capa más alta
                    delta{C}(k) = dsig(H{C-1}(k),B)*W{C}(l,k)*delta{C+1}(l); %ojo
                 end
             end
            
        end
        
        %actualizo pesos----------------------
        for C = 1:Nt-1 %avanza en capas
            
            for k = 1:N(C+1) %sig capa
                for l = 1:N(C) %capa actual
                    dW(k,l) = n*delta{C+1}(k)*V{C}(l);
                    W{C}(k,l) = dW(k,l) + W{C}(k,l);
                end
            end
            
        end
        O(indice(j)) = double(V{end});
    end 
    iter = iter + 1;
    ECM(iter) = (1/2) * (yd-O')'*(yd-O');
    
end
ECM = ECM(2:end);

%% Ploteo ECM

figure
plot(ECM)
title('Error cuadrático medio')
ylabel('ECM')
xlabel('Iteraciones')

%% Pruebo el algoritmo

A(1) = calc_salida(Nt, [0 0], W, N, B);
A(2) = calc_salida(Nt, [0 1], W, N, B);
A(3) = calc_salida(Nt, [1 0], W, N, B);
A(4) = calc_salida(Nt, [1 1], W, N, B);







