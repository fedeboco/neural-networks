close all;
clear all;

eta = 1;

x = [
    -1 -1 -1 -1;
    -1 -1 -1 1;
    -1 -1 1 -1;
    -1 -1 1 1;
    -1 1 -1 -1;
    -1 1 -1 1;
    -1 1 1 -1;
    -1 1 1 1;
    1 -1 -1 -1;
    1 -1 -1 1;
    1 -1 1 -1;
    1 -1 1 1;
    1 1 -1 -1;
    1 1 -1 1;
    1 1 1 -1;
    1 1 1 1
    ];

yd = [
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1;
    1
    ];

yd = fliplr(yd) .* -1; %comentar para and

nro_entradas = 4;
entrada = [x ones(2^nro_entradas, 1)];

W = randn(1, nro_entradas + 1);

aprendio = 0;
while ~aprendio
    indices = randperm(2^nro_entradas);
    for i = 1: length(entrada);
        y(indices(i)) = signo(W * entrada(indices(i), :)');
        delta_W = eta * entrada(indices(i), :) * (yd(indices(i)) - y(indices(i)));
        W = W + delta_W;
    end
    aprendio = isequal(yd, y');
    
end
