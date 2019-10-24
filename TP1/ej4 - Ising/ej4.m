%Federico Pérez Boco 96777 - 2do cuat. 2018

%Modelo de ISING, dipolos magnéticos.
%Resulta de una analogía entre redes neuronales y dipolos magnéticos, los
%cuales son dependientes con la temperatura.

close all;
clear all;

N=10; %dimension de matriz de dipolos
S_inicial = signo(randn(N)); %estados
S = S_inicial;

%Calculo Energia

E = calc_energia(S,N);
T0 = 5;
paso = 0.02;
k = 1;
ze = 1;

%Uso SN(I(i),J(i)) por cada par en i tiene la ubicación específica
%de un elemento permutado.
T=T0;
while(T > 0.01)
   for actualizaciones = 1:10 %veces que se actualiza cada dipolo
        [I,J] = ind2sub([N], randperm(N*N)); %indices aleatorios
        for i=1:N*N
                SN = S; %estados nuevos
                SN (I(i),J(i)) = -S (I(i),J(i));
                EN = calc_energia(SN,N); %energia nueva
                deltaE = EN - E;
                if(deltaE < 0) %bajó energía, acepto el cambio
                    S = SN;
                    E = EN;
                else
                    P = exp(-deltaE/(k*T));
                    if(rand < P) %acepto el cambio
                        S = SN;
                        E = EN;
                    end
                end      
        end
    end
    M(ze) = sum(sum(S));
    ze = ze + 1;
    T = T - paso
end
M = fliplr(M);
eje_T = paso:paso:T0;
figure
plot(eje_T, M)
xlabel('Temperatura [K]')
ylabel('Magnetización')

figure
subplot(1,2,1)
pcolor(S_inicial)
title('Estados iniciales')
subplot(1,2,2)
pcolor(S)
title('Magnetización final')