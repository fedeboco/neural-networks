close all;
clear all;

eta = 1;

x = [
    -1 -1;
    -1 1;
    1 -1;
    1 1
    ];

yd = [
    -1;
    -1;
    -1;
    1
    ];

yd = fliplr(yd) .* -1; %comentar para and

entrada = [x ones(4, 1)];

W = randn(1, 3);

aprendio = 0;
j = 0;
while ~aprendio
    indices = randperm(4);
    for i = 1: length(entrada);
        y(indices(i)) = signo(W * entrada(indices(i), :)');
        delta_W = eta * entrada(indices(i), :) * (yd(indices(i)) - y(indices(i)));
        W = W + delta_W;
    end
    
    aprendio = isequal(yd, y');
    
end

for i = 1: length(y)
   
    if y(i) == 1
        color(i, :) = [1 0 0];
    else
        color(i, :) = [0 0 0];
    end
    
end

x1 = [-2 2];
x2 = (-W(3) - W(1) .* x1) / W(2);

figure;
hold on;
scatter(x(:, 1), x(:, 2), [], color);
plot(x1, x2);
axis([-2 2 -2 2]);
