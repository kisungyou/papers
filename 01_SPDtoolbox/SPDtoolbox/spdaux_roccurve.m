function [vecTPR, vecFPR] = spdaux_roccurve(truelabel, pvalues)

nlength  = 1234;
myalphas = linspace(0.001,0.999, nlength);
vecTPR = zeros(1,nlength);
vecFPR = zeros(1,nlength);

tid0 = find(truelabel < 0.5);
tid1 = setdiff(1:length(truelabel), tid0);

for i=1:nlength
    alpha = myalphas(i);
    fakelabel = (pvalues < alpha)*1.0;
    fid1      = find(fakelabel > 0.5);
    
    vecTPR(i) = length(intersect(tid1,fid1))/length(tid1);
    vecFPR(i) = length(intersect(tid0,fid1))/length(tid0);
end
end