function salida = calc_salida(Nt, entrada, W, N, B)
    V{1} = [entrada 1];
    for C = 1:Nt-1 %avanza en Capas

        for k = 1:N(C+1) %calcula para cada neurona
            H{C}(k) = W{C}(k,:)*V{C}';
            V{C+1} = sig(H{C},B);
        end

    end
    salida = V{end};
end