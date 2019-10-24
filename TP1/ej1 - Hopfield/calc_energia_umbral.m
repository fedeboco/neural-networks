function E = calc_energia_umbral(s,W,umbral)
    %notar que el umbral podria ser un vector con escalares distintos para
    %cada neurona
    umbral = umbral.*ones(1,length(s));
    E = -1/2*s'*W*s + umbral*s;
end