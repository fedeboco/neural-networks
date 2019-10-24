function f = gdfact(x)
    f = 1;
    n = x;
    while(n>=1)
        f = f*n;
        n = n-1;
    end
end