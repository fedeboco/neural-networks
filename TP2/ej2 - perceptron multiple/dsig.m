function R = dsig(H, B)
    R = 2*B.*sig(H,B).*(1-sig(H,B));
end