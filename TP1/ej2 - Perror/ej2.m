%Federico Pérez Boco 96777 - 2do cuat. 2018

%Este algoritmo prueba que mientras mayor probabilidad de error se admite
%mayor es la cantidad de patrones que se pueden almacenar par un N dado.
%Debe verificarse:
%pmax/N = [0.105 0.138 0.185 0.37 0.61] para las probabilidades dadas.

N = 400; %neuronas por cada patrón

probabilidades = [0.001 0.0036 0.01 0.05 0.1];
patrones_max_pct = [0.01 0.01 0.09 0.2 0.3];

resultados = zeros(length(probabilidades), 2, 10);

for k = 1 : 10
    for i = 1 : length(probabilidades)
        m = 0; %número de patrones
        P_error = 0;
        while ((P_error < probabilidades(i)) && m <= N)
            
            %Genero m patrones de N neuronas
            m = m + 1;
            P = sign(randn(N, m)); %Patrón: vector aleatorio
            W = (P * P') - m * eye(N); %actualizo matriz estados
            
            %Checkeo error
            errores = 0;
            for j = 1 : m
                resultado = signo(W * P(:, j));
                errores = errores + sum(resultado ~= P(:, j));
            end
            P_error = errores / (N * m);
            
        end
        resultados(i, :, k) = [P_error (m / N)];
    end
end

media = zeros(length(probabilidades), 2);

for i = 1 : length(probabilidades)
    media(i, :) = [mean(resultados(i, 1, :)) mean(resultados(i, 2, :))];
end

media
