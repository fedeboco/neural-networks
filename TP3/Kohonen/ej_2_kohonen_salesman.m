close all
clear all

%% Configuración
N = 2000; %nro de puntos por dimensión
CIUDADES = 200;
n = 0.8; %cte aprendizaje
s_iter = 300; %iteraciones con el mismo sigma
s0 = 20; %sigma inicial
sf = 5; %sigma final
iMAX = N/s_iter; %cant. máxima de iteraciones

%% Genero distribución ciudades

x=10*randn(1,CIUDADES);
y=10*randn(1,CIUDADES);

scatter(x,y,'.')

%% Entreno matriz de pesos sinápticos

iMAX = 1000*CIUDADES;
ang=2*pi*linspace(0,1,CIUDADES*2+1);
wx=cos(ang);
wy=sin(ang);
s = linspace(CIUDADES, 0.1, iMAX);	
iteracion = 0;

%ciudad inicial
wx(1)=x(1);
wy(1)=y(1);

indice = 2:(CIUDADES*2+1);

%para cada muestra
while (iteracion < iMAX)
    for m = 1:CIUDADES
        
        %cierro anillo
        wx(CIUDADES*2+2)=wx(1);
        wy(CIUDADES*2+2)=wy(1);

        %ajusto sigma
        sigma = s(iteracion+1);

        %busco ganadora
        d = sqrt(((x(m)-wx).^2) + ((y(m)-wy).^2)); % distancia
        [Igan, Jgan] = find(d==min(min(d)));
        G = [Igan(1) Jgan(1)]; %ganadora, (1) por si hay más de 1 mínimo

        %actualizo vecindad
        for i = 2 : CIUDADES*2+1
                dist = sqrt((i - G(1))^2);
                V = exp(-1.*((dist).^2)./(2*(sigma^2))); %vecindad
                dwx(i) = n*V.*(x(m)-wx(i));
                dwy(i) = n*V.*(y (m)-wy(i));
        end

        W{m,1} = wx; %para plotear
        W{m,2} = wy; %para plotear

        %actualizo matriz de pesos
        wx(2 : CIUDADES*2+1) = wx(2 : CIUDADES*2+1) + dwx(2 : CIUDADES*2+1);
        wy(2 : CIUDADES*2+1) = wy(2 : CIUDADES*2+1) + dwy(2 : CIUDADES*2+1);
    iteracion = iteracion + 1;
    end   
end
%% Plot

figure

plot(x,y,'ok')
hold on
plot(wx,wy)