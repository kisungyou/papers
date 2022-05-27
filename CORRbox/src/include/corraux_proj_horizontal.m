function horiV = corraux_proj_horizontal(S, V)

% S in SPD(n)
% V in T_S SPD(n)

horiV = V - corraux_proj_vertical(S,V);


end