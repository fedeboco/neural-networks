%Federico Pérez Boco 96777 - 2do cuat. 2018

%%Cargo imágenes binarias (Blanco y Negro) y normalizo a -1 y 1

i1 = cargar_patron('i1.bmp');
i2 = cargar_patron('i2.bmp');
i3 = cargar_patron('i3.bmp');
ir1 = cargar_patron('ir1.bmp');
ir2 = cargar_patron('ir2.bmp');
ir3 = cargar_patron('ir3.bmp');

%% Aplico algoritmo de Hopfield

P = [i1, i2, i3]; %Matriz de patrones
N = 3; %Numero de patrones
L = length(i1); %Numero de neuronas
W = (P * P') - N * eye(L); %Matris de pesos sincrónica

%% Verifico que las imagenes originales son de energía mínima local

r = sign(W * i1);
d1 = sum(r~=i1);

r = sign(W * i2);
d2 = sum(r~=i2);

r = sign(W * i3);
d3 = sum(r~=i3);

if(d1~=0 || d2~=0 || d3~=0)
    display('No es de energía mínima.')
end    

%% Filtro imagenes con ruido

%patron_act = actualizar_asincronico(mat_pesos, patron_ruidoso, iteraciones, plotear progreso (bool));
r = actualizar_asincronico(W, ir1, 1, 1);
d1 = sum(r~=i1)
reconst1 = (vec2mat(r,sqrt(L)))';

r = actualizar_asincronico(W, ir2, 1, 1);
d2 = sum(r~=i2)
reconst2 = (vec2mat(r,sqrt(L)))'; 

r = actualizar_asincronico(W, ir3, 1, 1);
d3 = sum(r~=i3)
reconst3 = (vec2mat(r,sqrt(L)))';  

imwrite(reconst1,'reconst1.bmp')
imwrite(reconst2,'reconst2.bmp')
imwrite(reconst3,'reconst3.bmp')

%% Pruebo con imagenes invertidas
i1_inv = -i1;
i2_inv = -i2;
i3_inv = -i3;
r = actualizar_asincronico(W, i1_inv, 1, 1);
reconst1_inv = (vec2mat(r,sqrt(L)))';
r = actualizar_asincronico(W, i2_inv, 1, 1);
reconst2_inv = (vec2mat(r,sqrt(L)))';
r = actualizar_asincronico(W, i3_inv, 1, 1);
reconst3_inv = (vec2mat(r,sqrt(L)))';

%% Pruebo con mezclas
i1(i1 == -1) = 0;
i2(i2 == -1) = 0;
i3(i3 == -1) = 0;

i1_i3 = double(i1 | i3);
i2_i3 = double(i2 & i3);
i3_i1 = double(i1 & i3);

i1_i3(i1_i3 == 0) = -1;
i2_i3(i2_i3 == 0) = -1;
i3_i1(i3_i1 == 0) = -1;

r = actualizar_asincronico(W, i1_i3, 1, 1);
reconst1_mezc = (vec2mat(r,sqrt(L)))';
r = actualizar_asincronico(W, i2_i3, 1, 1);
reconst2_mezc = (vec2mat(r,sqrt(L)))';
r = actualizar_asincronico(W, i3_i1, 1, 1);
reconst3_mezc = (vec2mat(r,sqrt(L)))';

%% Ploteo

figure(1)
subplot(3,3,1)
imshow(vec2mat(i1, sqrt(L))');
title('Original 1')
subplot(3,3,2)
imshow(vec2mat(i2, sqrt(L))')
title('Original 2')
subplot(3,3,3)
imshow(vec2mat(i3, sqrt(L))')
title('Original 3')

subplot(3,3,4)
imshow('ir1.bmp')
title('Ruidosa 1')
subplot(3,3,5)
imshow('ir2.bmp')
title('Ruidosa 2')
subplot(3,3,6)
imshow('ir3.bmp')
title('Ruidosa 3')

subplot(3,3,7)
imshow(reconst1)
title('Reconstruida 1')
subplot(3,3,8)
imshow(reconst2)
title('Reconstruida 2')
subplot(3,3,9)
imshow(reconst3)
title('Reconstruida 3')

%% Ploteo invertidas

figure(2)
subplot(3,3,1)
imshow(vec2mat(i1, sqrt(L))');
title('Original 1')
subplot(3,3,2)
imshow(vec2mat(i2, sqrt(L))')
title('Original 2')
subplot(3,3,3)
imshow(vec2mat(i3, sqrt(L))')
title('Original 3')

subplot(3,3,4)
imshow(vec2mat(i1_inv, sqrt(L))')
title('Invertida 1')
subplot(3,3,5)
imshow(vec2mat(i2_inv, sqrt(L))')
title('Invertida 2')
subplot(3,3,6)
imshow(vec2mat(i3_inv, sqrt(L))')
title('Invertida 3')

subplot(3,3,7)
imshow(reconst1_inv)
title('Reconstruida 1')
subplot(3,3,8)
imshow(reconst2_inv)
title('Reconstruida 2')
subplot(3,3,9)
imshow(reconst3_inv)
title('Reconstruida 3')

%% Ploteo mezclas

figure(3)
subplot(3,3,1)
imshow(vec2mat(i1, sqrt(L))');
title('Original 1')
subplot(3,3,2)
imshow(vec2mat(i2, sqrt(L))')
title('Original 2')
subplot(3,3,3)
imshow(vec2mat(i3, sqrt(L))')
title('Original 3')

subplot(3,3,4)
imshow(vec2mat(i1_i3, sqrt(L))')
title('Mezcla AND 1 y 3')
subplot(3,3,5)
imshow(vec2mat(i2_i3, sqrt(L))')
title('Mezcla OR 2 y 3')
subplot(3,3,6)
imshow(vec2mat(i3_i1, sqrt(L))')
title('Mezcla OR 1 y 3')

subplot(3,3,7)
imshow(reconst1_mezc)
title('Reconstruida 1')
subplot(3,3,8)
imshow(reconst2_mezc)
title('Reconstruida 2')
subplot(3,3,9)
imshow(reconst3_mezc)
title('Reconstruida 3')