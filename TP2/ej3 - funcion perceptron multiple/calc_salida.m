function salida = calc_salida(Nt, entrada, W, N, B)
    
    V{1} = entrada;
    for C = 1:Nt-1 %avanza en capas C=CAPA
        for k = 1:N(C+1)-1 %calcula para cada neurona excepto umbral

            H{C}(k) = W{C}(k,:)*[V{C} 1]'; %suma para neurona k-esima, 1 por umbral
        end
        V{C+1} = sig2(H{C},B); %sigmoide en todas las neuronas de la capa
    end
    salida = V{end};
end