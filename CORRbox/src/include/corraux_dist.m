function dval = corraux_dist(A,B)
    D = corraux_findD(A,B);
    C = D*B*D;
    dval = corraux_spddist(A,C);
end