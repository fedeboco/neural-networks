function patron = cargar_patron(archivo)
    patron = imread(archivo); 
    patron = double(patron(:, :, 1)./255);
    patron(patron == 0) = -1;
    patron = patron(:);

end