
%Federico Pérez Boco 96777 - Redes Neuronales

clear all
close all
%Págs 116 y 120

%% Entrada

%patrón de entrada
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
            yd = [yd; (sin(i) + cos(j) + k)./3]; %salida deseada; ./3 por tanh
        end
    end
end

 
%% Configuración
 
%Número de entradas
Li = length(p(1,:));

%Número de patrones
Lp = length(p(:,1));

%Número de salidas
Lo = length(yd(1,:));

%Elementos de capas i-ésima (1 capa de entradas, k capas de neuronas, 1 capa de salidas)
N = [Li 5 2 Lo] + 1; %el +1 es por los umbrales en cada capa
Nt = length(N); %Número de elementos/Cs

%Beta función sigmoide
B = 1/2; %tipico 1/2 o 1

%Constante de aprendizaje nu
n = 0.1; %entre 0 y 1

%Inicializo pesos con valores pequeños
for i = 1:Nt-1
    W{i} = rand(N(i+1)-1,N(i))./100;  
end

%Umbral de ECM
U = 1;

%Iteraciones máximas
iter_MAX = 1000;

%% Procesamiento

ECM(1) = inf;
iter = 1;

while(ECM(iter) > U || iter < iter_MAX)
    
    indice = randperm(length(p(:,1))); %indice aleatorio para patrones
    ECM_parcial = 0;
    for j = 1:Lp %permuta entre patrones

        %forward-----------------------------
        V{1} = [p(indice(j),:)];
        for C = 1:Nt-1 %avanza en capas C=CAPA
            for k = 1:N(C+1)-1 %calcula para cada neurona excepto umbral

                H{C}(k) = W{C}(k,:)*[V{C} 1]'; %suma para neurona k-esima, 1 por umbral
            end
            V{C+1} = sig2(H{C},B); %sigmoide en todas las neuronas de la capa
        end
        
        %backward propagation----------------
        delta{Nt} = dsig2(H{end},B).*(yd(indice(j),:)-V{end});
        for C = Nt-1:-1:2 %para cada capa excepto la última

             for k = 1:N(C)-1 %capa actual
                 for l = 1:N(C+1)-1 %capa más alta
                    delta{C}(k) = dsig2(H{C-1}(k),B)*W{C}(l,k)*delta{C+1}(l);
                 end
             end
            
        end
        
        %actualizo pesos----------------------
        for C = 1:Nt-1 %avanza en capas            
            
            dW{C} = n*delta{C+1}'*[V{C} 1];
            W{C} = dW{C} + W{C};
                        
        end
        
        %actualizo ECM------------------------
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
title('Error cuadrático medio')
ylabel('ECM')
xlabel('Iteraciones')

%% Pruebo el algoritmo

clear y_test xtest ptest g

puntos_test = 100;
x_test = 1:puntos_test;

%prueba 1: sin(x)+cos(0)+0 ------------------------------

xtest = linspace(0,2*pi,puntos_test);

ptest = [xtest; zeros(1,puntos_test); zeros(1,puntos_test)]';
for i=1:length(ptest)
    g(i,:) = calc_salida(Nt, ptest(i,1:end), W, N, B);
end

figure
subplot(1,3,1)
plot(xtest,  g)
hold on
plot(xtest, sin(xtest)./3,'--')
stem(x, sin(x)./3)
legend('aproximada', 'real', 'muestras')

%prueba 2: sin(0)+cos(x)+0 ------------------------------

xtest = linspace(0,2*pi,puntos_test);
ptest = [zeros(1,puntos_test); xtest; zeros(1,puntos_test)]';
for i=1:length(ptest)
    g(i,:) = calc_salida(Nt, ptest(i,1:end), W, N, B);
end


subplot(1,3,2)
plot(xtest,  g)
hold on
plot(xtest, cos(xtest)./3,'--')
stem(y, cos(x)./3)
legend('aproximada', 'real', 'muestras')

%prueba 3: sin(0)+cos(0)+x ------------------------------

xtest = linspace(-1,1,puntos_test);
ptest = [zeros(1,puntos_test); zeros(1,puntos_test); xtest]';
for i=1:length(ptest)
    g(i,:) = calc_salida(Nt, ptest(i,1:end), W, N, B);
end


subplot(1,3,3)
plot(xtest,  g)
hold on
plot(xtest, xtest./3,'--')
stem(z, z./3)
legend('aproximada', 'real', 'muestras')