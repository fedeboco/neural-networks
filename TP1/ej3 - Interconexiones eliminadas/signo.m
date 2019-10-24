function s = signo(x)
    s = sign(x);
    s(s == 0) = 1;
end
