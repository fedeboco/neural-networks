function patron = actualizar_asincronico(W, patron, iteraciones, plot_figure)
    L = length(patron);
    v = randperm(L); %genero indices aleatorios
    umbral = 0;
    for z =1:iteraciones
        for i=1:L
            suma = W(v(i),:)*patron;
            if(suma >= umbral)
                patron(v(i)) = 1;
            else
                patron(v(i)) = -1;
            end
            figure_plotted = 0;
            if(plot_figure)
                figure(50)
                imshow( (vec2mat(patron,sqrt(L)))', 'InitialMagnification', 400 )
                figure_plotted = 1;
            end
        end
    end
    if(figure_plotted)
        close(50)
    end
end