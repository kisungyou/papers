%% AUXILIARY FUNCTIONS
function output = corraux_test2bg_statistic(DX,DY,DXY)
m = size(DXY,1);
n = size(DXY,2);

muff = sum(sum(triu(DX,1)))/(m*(m-1)/2);
mufg = sum(sum(DXY))/(m*n);
mugg = sum(sum(triu(DY,1)))/(n*(n-1)/2);

vec1 = [muff,mufg];
vec2 = [mufg,mugg];
output = sum((vec1-vec2).^2);
end