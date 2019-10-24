%Federico Pérez Boco 96777 - 2do cuat. 2018
close all
clear all

N = 100; %nro de neuronas
m = 10; %nro de patrones

P_errores = [];

P = signo(randn(N, m));
W = (P * P') - m * eye(N);

for k = 0 : 100 % cantidad de conexiones eliminadas (en %)
    A = (rand(N) > (0.01 * k)); %asigno peso 0 (elimina conexión) aleatorio
    W_nueva = W .* A;
    errores = 0;
    for j = 1 : m
        resultado = signo(W_nueva * P(:, j));
        errores = errores + sum(resultado ~= P(:, j));
    end
    P_error = errores / (N * m);
    P_errores = [P_errores, P_error];
end

figure(1);
plot(0:1:100, P_errores);
ylabel('p_{error}')
xlabel('Conexiones eliminadas [%]')