function e = gdexp(x)
    L = length(x);
    e = ones(1,L);
    for t = 1:L
        for i=1:50
            z = 1;
            for j=1:i
                z = z*x(t);
            end
            f = gdfact(i);
            e(t) = e(t) + z/f;
        end
    end
end