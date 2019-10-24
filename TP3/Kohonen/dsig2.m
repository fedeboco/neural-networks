function R = dsig2(H, B)
    R = B.*(1-sig2(H,B).^2);
end