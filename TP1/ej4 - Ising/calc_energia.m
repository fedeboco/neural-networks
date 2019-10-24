function E = calc_energia(S,N)

    %rodeo con ceros para poder sumar vecinos sin problemas
    S = [zeros(length(S),1) S zeros(length(S),1)];
    S = [zeros(1,length(S)); S; zeros(1,length(S))];

    E = 0;
    J = 1;
    for i = 2:N+1
        for j = 2:N+1
            E = E + S(i,j)*(S(i-1,j)+S(i+1,j)+S(i,j-1)+S(i,j+1));
        end
    end
    E = J*0.5*E;
end