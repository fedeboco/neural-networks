clear all
close all

%% Número de ciudades
CIUDADES=200;

%% Anillo de neuronas
indice = 2:(CIUDADES*2+1); 
phi = 2*pi*linspace(0, 1, CIUDADES*2+1);
wx = cos(phi);
wy = sin(phi);

%ploteo anillo
figure
plot(wx,wy,'b')
hold on
plot(wx,wy,'.r')

%% Ciudades aleatorias
x = 10*randn(1,CIUDADES);
y = 10*randn(1,CIUDADES);

%ciudad inicial
wx(1)=x(1);
wy(1)=y(1);

%% Configuración
n = 0.8; %cte de aprendizaje
IMAX = 100000; %MAX ITERACIONES
s = linspace(CIUDADES, 0.1, IMAX); %ancho de vecindad

%% Procesamiento
for i=1:IMAX
    
    j = randi(CIUDADES); %ciudad al azar
    
    wx(CIUDADES*2+2) = wx(1); % cierro anillo
    wy(CIUDADES*2+2) = wy(1); % cierro anillo
    
    d = sqrt( (wx-x(j)).^2 + (wy-y(j)).^2 ); %distancias
	[~, G] = min(d); %indice ganador de las distancias
    
    sigma = s(i); %de la vecindad para iteración i-esima
    V = exp(-(abs(indice-G).^2)/(2*(sigma^2))); %vecindad
    
    %actualizo pesos
	wx(2:CIUDADES*2+1) = wx(2:CIUDADES*2+1) + n*V.*(x(j) - wx(2:CIUDADES*2+1));
	wy(2:CIUDADES*2+1) = wy(2:CIUDADES*2+1) + n*V.*(y(j) - wy(2:CIUDADES*2+1));

end

%% Ploteo final

figure
plot(wx,wy,'r')
hold on
plot(x,y,'ok')
hold off
title('Salesman problem')