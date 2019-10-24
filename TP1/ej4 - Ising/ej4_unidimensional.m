%Federico Pérez Boco 96777 - 2do cuat. 2018

%Modelo de ISING, dipolos magnéticos.
%Resulta de una analogía entre redes neuronales y dipolos magnéticos, los
%cuales son dependientes con la temperatura.

close all;
clear all;

N=100; %dimension de matriz de dipolos
S_inicial = signo(randn(N,1)); %estados
S = S_inicial;

%Calculo Energia

E = 0;
J = 1;
S = [0; S; 0];
for i = 2:N+1
        E = E + S(i)*(S(i-1)+S(i+1));
end
E = J*0.5*E;
S = S(2:N+1);

T0 = 5;
paso = 0.02;
k = 1;
ze = 1;

%Uso SN(I(i),J(i)) por cada par en i tiene la ubicación específica
%de un elemento permutado.
T=T0;
while(T > 0.01)
   for actualizaciones = 1:10 %veces que se actualiza cada dipolo
        I = randperm(N); %indices aleatorios
        for i=1:N
                SN = S; %estados nuevos
                SN (I(i)) = -SN (I(i));
                SN = [0; SN; 0];
                EN=0;
                for j = 2:N+1
                    EN = EN + SN(j)*(SN(j-1)+SN(j+1));
                end
                EN = J*0.5*EN; %energia nueva
                deltaE = EN - E;
                SN = SN(2:N+1);
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
    M(ze) = sum(S);
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
pcolor(vec2mat(S_inicial, sqrt(N)))
title('Estados iniciales')
subplot(1,2,2)
pcolor(vec2mat(S, sqrt(N)))
title('Magnetización final')