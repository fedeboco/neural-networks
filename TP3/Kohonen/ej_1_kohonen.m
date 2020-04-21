%Federico Pérez Boco 96777 - Redes Neuronales

close all
clear all

%% Configuración
N = 2000; %nro de puntos por dimensión
Nn = 20; %nro de neuronas
n = 0.1; %cte aprendizaje
s_iter = 300; %iteraciones con el mismo sigma
s0 = 20; %sigma inicial
sf = 5; %sigma final
iMAX = N/s_iter; %cant. máxima de iteraciones

%% Genero distribución uniforme en círculo unitario
%http://mathworld.wolfram.com/DiskPointPicking.html

%circulo
% r = rand(1,N); %uniforme [0,1]
% t = rand(1,N)*2pi; %uniforme [0,2pi)
% x = sqrt(r).*cos(t); %a circ unitario
% y = sqrt(r).*sin(t); %a circ unitario

%triangulo
j = 0;
x = zeros(1,N);
while( j < N)
    auxx = rand;
    auxy = rand;
    if(auxy < 2*auxx && auxy < -2*auxx+1)
        j = j + 1;
        x(j) = auxx;
        y(j) = auxy;
    end
end

figure
scatter(x,y,'.g')

%% Entreno matriz de pesos sinápticos

wx = 2*rand(Nn,Nn)-1;
wy = 2*rand(Nn,Nn)-1;
s = logspace(log10(sf), log10(s0), iMAX + 1);	
s = fliplr(s);
t = 0; %iteración
is = 0;
ix = 1:Nn;
iy = 1:Nn;


%para cada muestra
for m = 1:N

    %ajusto sigma
    sigma = s(ceil(m / s_iter + 0.001));

    %busco ganadora
    d = sqrt(((abs(x(m)-wx)).^2) + ((abs(y(m)-wy)).^2)); % distancia
    [Igan, Jgan] = find(d==min(min(d)));
    G = [Igan(1) Jgan(1)]; %ganadora, (1) por si hay más de 1 mínimo
    
    %actualizo vecindad
    for i = 1 : Nn
        for j = 1 : Nn
            dist = sqrt((i - G(1))^2 + (j - G(2))^2);
            V = exp(-1.*((dist).^2)./(2*(sigma^2))); %vecindad
            dwx(i, j) = n*V.*(x(m)-wx(i,j));
            dwy(i, j) = n*V.*(y(m)-wy(i,j));
        end
    end
    
    W{m,1} = wx; %para plotear
    W{m,2} = wy; %para plotear
    
    %actualizo matriz de pesos
    wx = wx + dwx;
    wy = wy + dwy;
    
    %aumento contadores
    t = t + 1; %contador de iteracion
    is = is + 1; %contador de iteraciones con mismo sigma
    m = m + 1; %contador de muestra

end
  

%% Plot

figure

subplot(2,2,1)
hold on;
plot(W{1,1},  W{1,2}, 'r', 'LineWidth', 1)
plot(W{1,1}', W{1,2}','r', 'LineWidth', 1)
plot(W{1,1},  W{1,2}, 'ok');
axis([-1 1 -1 1])
title('Primera iteración');
xlabel('w_x')
ylabel('w_y')

subplot(2,2,2)
hold on;
plot(W{floor(end/50),1},  W{floor(end/50),2}, 'r', 'LineWidth', 1)
plot(W{floor(end/50),1}', W{floor(end/50),2}','r', 'LineWidth', 1)
plot(W{floor(end/50),1},  W{floor(end/50),2}, 'ok');

title('2%');
xlabel('w_x')
ylabel('w_y')

subplot(2,2,3)
hold on;
plot(W{floor(end/30),1},  W{floor(end/30),2}, 'r', 'LineWidth', 1)
plot(W{floor(end/30),1}', W{floor(end/30),2}','r', 'LineWidth', 1)
plot(W{floor(end/30),1},  W{floor(end/30),2}, 'ok');

title('3.3%');
xlabel('w_x')
ylabel('w_y')

subplot(2,2,4)
hold on;
plot(wx,  wy, 'r', 'LineWidth', 1)
plot(wx', wy','r', 'LineWidth', 1)
plot(wx,  wy, 'ok');
axis([-1 1 -1 1])
title('Última iteración');
xlabel('w_x')
ylabel('w_y')
