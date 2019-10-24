function V = sig(H, B)
    V = 1 ./ (1 + exp(-2*B*H));
end