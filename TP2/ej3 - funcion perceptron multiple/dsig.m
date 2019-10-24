function R = dsig(H, B)
    alfa(:,:) = sig(H,B);
    R = 2*B.*alfa.*(1-alfa);
end