%Federico Pérez Boco 96777 - Redes Neuronales - FIUBA

%Este algoritmo va variando la matriz de pesos simulando las variaciones
%geneticas en una poblacion. Se busca resolver el problema de obtener la
%salida deseada yd encontrando las mejores matrices W y w. Para ello se
%generan 100 pares de matrices aleatorias W y w (100 individuos) y se
%escogen las mejores 10. A las otras 90 primero se les mezcla los genes
%(valores de W y w) tomando algunos los valores de W y w de un individuo y
%poniéndoselos a otro y viceversa. Luego se les aplica un factor de mutación
%a los 100 individuos y se calcula la funcion de fitness que es un
%indicador de que tan bien se adecua el individuo a la salida. Mientras mas
%alta mejor, es la inversa del ECM. Luego calculo las probabilidades de
%reproducirse de cada individuo y obtengo la siguiente generacion en base a
%ellas. Finalmente calculo el fitness de todos y me quedo con el mejor, o
%sea las W y w con menor ECM.

close all;
clear all;

%% Configuracion
%Patrones
x = [ -1 -1; -1  1; 1 -1; 1  1]; %entradas
x = [x, ones(length(x), 1)]; %agrego umbral
yd = [ -1 1 1 -1 ]; %salida deseada 

%Propiedades de red
Li = size(x,2) - 1; %cantidad de entradas (sin umbral)
N1 = 2; %neuronas capa 1
N2 = 1; %neuronas capa 2
beta = 1; % beta
Ni = 100; %cantidad de individuos
Dm = 20; %desvio de mutacion
Pm = 0.01; %probabilidad de mutacion
U = 0.01; %umbral de ECM
P_reproduccion_inicial = 0.001; %Probabilidad de reproducción inical

nC1 = N1*(Li + 1); %numero de elementos de W en primer capa
nC2 = N2*(N1 + 1) + nC1; %numero de elementos de W y w en ambas capas

%% Procesamiento

%Generacion inicial
individuos = [];
for i = 1 : Ni
    
    W = randn(N1, Li + 1); %capa 1
    w = randn(N2, N1 + 1); %capa 2
    individuos = [individuos [W(:); w(:)]]; %paso a vector y concateno
    
end

ECM(1) = inf;
t = 0;
while(ECM > U)
    
    t = t + 1;
    
    %Calculo fitness
    F = zeros(1, Ni); %fitness inicial
    for i = 1 : Ni
        
        W = vec2mat(individuos(1 : nC1, i), Li + 1); %primer capa
        w = vec2mat(individuos(nC1 + 1 : nC2, i), N1 + 1); %segunda capa
        F(i) = fitness(individuos(:,1), W, w, x, yd, Li, beta);
        
    end
    [~, I] = sort(F, 'descend');
    individuos = individuos(:, I); %ordeno decreciente
    
    %Mezclo genes del 90% de los individuos con peor fitness
    K = round(Ni/10);
%     nueva_gen = zeros(size(individuos,1), size(individuos,2));
    nueva_gen = individuos(:, 1:K); %hijos con mejor fitness
    for i = 1 : Ni-K
        
        padre = randi(Ni); %padre
        madre = randi(Ni); %madre
        corte = randi(nC2); %cantidad de genes al azar que conservo   
        genesP = individuos(1 : corte, padre); %genes del padre
        genesM = individuos(corte + 1 : nC2, madre); %genes de la madre        
        nuevo_individuo = [genesP; genesM]; %hijo
        nueva_gen = [nueva_gen nuevo_individuo]; %nueva generacion
        
    end
    individuos = nueva_gen;
    
    %Muto la nueva generacion completa
    for i = 1 : nC2
        for j = 1 : Ni
            p = rand;
            if(p < Pm)
                individuos(i, j) = individuos(i, j) + randn*Dm;
            end
        end
    end
    
    %Nuevo fitness
    F = zeros(1, Ni);
    for i = 1 : Ni
        W = vec2mat(individuos(1 : nC1, i), Li + 1);
        w = vec2mat(individuos(nC1 + 1 : nC2), N1 + 1);
        F(i) = fitness(individuos, W, w, x, yd, Li, beta);
    end
    
    %Ordeno individuos y fitness segun fitness
    [F, I] = sort(F, 'descend');
    individuos = individuos(:, I);
    
    %Reproducción
    F_total = sum(F); %fitness total
    P = F ./ F_total; %probabilidad de reproducirse de cada indiv.
    nueva_gen = zeros(size(individuos,1), size(individuos,2));
    for i = 1 : Ni
        
        j = 1; %indice de individuo de gen anterior que pasa a ser de gen nueva
        p = rand; %probabilidad aleatoria
        P_acum = P_reproduccion_inicial; %probabilidad acumulada, empieza muy chica
        
        %Creo nueva generación
        while(j <= Ni)
            if(p < P_acum)
                nueva_gen(:,i) = individuos(:, j); %busco nuevo individuo i-esimo
                break;
            else
                if(j == 100)
                    j = 0;
                end
                j = j + 1;
                P_acum = P_acum + P(j); %aumento prob de reproducirse
            end
        end
        
    end
    individuos = nueva_gen;
    
    %Calculo fitness
    F = zeros(1, Ni);
    for i = 1 : Ni
        W = vec2mat(individuos(1 : nC1, i), Li + 1);
        w = vec2mat(individuos(nC1 + 1 : nC2, i), N1 + 1);
        F(i) = fitness(individuos, W, w, x, yd, Li, beta);
    end
    
    %Ordeno individuos y fitness segun fitness
    [F, I] = sort(F, 'descend');
    individuos = individuos(:, I);
    
    %El mejor elemento
    ECM(t) = 1 / fitness(individuos(:, 1), W, w, x, yd, Li, beta);
    ECM(t)
end

%% Pruebo salida

W = vec2mat(individuos(1 : nC1, 1), Li + 1); %1 por ser el mejor individuo
w = vec2mat(individuos(nC1 + 1 : nC2, 1), N1 + 1);

y = zeros(2^Li, 1);
for i = 1:2^Li
    h1 = W * x(i, :)';
    v = tanh(beta * h1);
    h2 = w * [v; 1];
    y(i) = tanh(beta * h2);
end
     
ECM
y

%% Ploteo

x1 = [];
x2 = [];
for i = 1:2^Li
   if(yd(i) == 1)
       x1 = [x1; x(i,:)];
   else
       x2 = [x2; x(i,:)];
   end
end

figure;
plot(ECM)

figure

hold on
plot(x1(:,2), x2(:,2), 'xb');
plot(x1(:,1), x2(:,1), 'xr');
axis([-1.5, 1.5, -1.5, 1.5]);

discx = -1.5:0.01:1.5;
discy = - W(1,1) / W(1,2) * discx - W(1,3) / W(1,2);
plot(discy, discx, 'k--');
discy = - W(2,1) / W(2,2) * discx - W(2,3) / W(2,2);
plot(discy, discx, 'k--');

legend('Salida = 1', 'Salida = -1', 'Discriminante','location','best');