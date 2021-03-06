
%Federico P�rez Boco 96777 - Redes Neuronales

clear all
close all
%P�gs 116 y 120

%% Entrada

%XOR

%patr�n de entrada
puntos = 10;
x = (0:puntos-1)*2*pi/(puntos-1);
y = x;
z = (0:puntos-1)*2/(puntos-1) - 1;

pdiag = [x; y; z]'; %lo uso para testing
ydiag = [sin(x) + cos(y) + z];

p = [];
yd = [];

for i = x
    for j = y
        for k = z
            p = [p; i j k];
            yd = [yd; (sin(i) + cos(j) + k)./3]; %salida deseada
        end
    end
end

 
%% Configuraci�n
 
%N�mero de entradas
Li = length(p(1,:));

%N�mero de patrones
Lp = length(p(:,1));

%N�mero de salidas
Lo = length(yd(1,:));

%Elementos de capas i-�sima (1 capa de entradas, k capas de neuronas, 1 capa de salidas)
N = [Li+1 15+1 Lo+1];
Nt = length(N); %N�mero de elementos/Cs

%Beta funci�n sigmoide
B = 1/2; %tipico 1/2 o 1

%Constante de aprendizaje nu
n = 0.1; %entre 0 y 1

%Inicializo pesos con valores peque�os
for i = 1:Nt-1
    W{i} = rand(N(i+1),N(i))./100;  
end

%Umbral de ECM
U = 0.01;

%% Procesamiento

O = zeros(1,4);
ECM(1) = inf;
iter = 1;

while(ECM(iter) > U)
    
    indice = randperm(length(p(:,1))); %indice aleatorio para patrones
    ECM_parcial = 0;
    for j = 1:Lp %permuta entre patrones

        %forward-----------------------------
        V{1} = [p(indice(j),:)];
        for C = 1:Nt-1 %avanza en capas C=CAPA
            
            for k = 1:N(C+1) %calcula para cada neurona
                H{C}(k) = W{C}(k,:)*[V{C} 1]'; %suma para neurona k-esima, 1 por umbral
            end
            V{C+1} = sig2(H{C},B); %sigmoide en todas las neuronas de la capa
            V{C+1} = V{C+1}(1:end-1); %elimino la ultima salida porque corresponde al umbral
            H{C} = H{C}(1:end-1);
        end
        
        %backward----------------------------
        delta{Nt} = dsig2(H{end},B).*(yd(indice(j),:)-V{end});
        for C = Nt-1:-1:2 %para cada capa excepto la �ltima

             for k = 1:N(C)-1 %capa actual
                 for l = 1:N(C+1)-1 %capa m�s alta
                    delta{C}(k) = dsig2(H{C-1}(k),B)*W{C}(l,k)*delta{C+1}(l);
                 end
             end
            
        end
        
        %actualizo pesos----------------------
        for C = 1:Nt-1 %avanza en capas
            
            for k = 1:N(C+1)-1 %neurona de la sig capa
                for l = 1:N(C)-1 %neurona de la capa actual
                    dW(k,l) = n*delta{C+1}(k)*V{C}(l);
                    W{C}(k,l) = dW(k,l) + W{C}(k,l); %!!!!!!!!!!!!!!!!!!!!! dW 15x15!
                end
            end
            
        end
        
        ECM_parcial = ECM_parcial + sum((yd(indice(j),:)-[V{end}])*(yd(indice(j),:)-[V{end}])');  
    end 
    iter = iter + 1;
    ECM(iter) = ECM_parcial;
    ECM(iter)

    
end
ECM = ECM(2:end);

%% Ploteo ECM

figure
plot(ECM)
title('Error cuadr�tico medio')
ylabel('ECM')
xlabel('Iteraciones')

%% Pruebo el algoritmo

close all

for i=1:length(pdiag)
    g(i,:) = calc_salida(Nt, pdiag(i,1:end), W, N, B);
end

g

figure
plot(ydiag)
hold on
plot(g(:,1))




