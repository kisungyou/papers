function C3 = corraux_interpolate(C1, C2, w)

if ((w <= 0)||(w>=1))
    error(" interpolate : a weight 'w' should be in (0,1).");
end

%   : it computes C3 = w*C1 + (1-w)*C2
C3 = corraux_exp(C1, corraux_log(C1, C2), 1-w);

end