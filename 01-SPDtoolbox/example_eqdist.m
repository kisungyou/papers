%% EXAMPLE : test equality of distributions : empirical Type 1 error

%   preliminary
clear; close; clc; addpath(genpath(pwd));

%   example setup
ntest = 100;
p     = 5;              % dimension of covariance matrices
recx  = zeros(1,ntest); % vector to record p-values via extrinsic distance
reci  = zeros(1,ntest); %                               intrinsic distance

n1 = 10; % number of objects for class 1
n2 = 8;  %                       class 2

%   iterate !
for i=1:ntest
    dat1 = zeros(p,p,n1);
    dat2 = zeros(p,p,n2); 
    
    for m=1:n1
        dat1(:,:,m) = cov(randn(100,p));
    end
    for n=1:n2
        dat2(:,:,n) = cov(randn(100,p));
    end
    
    spdobj1 = spd_initialize(dat1); % don't forget to initialize the data !
    spdobj2 = spd_initialize(dat2); 
    
    recx(i) = spd_eqtest(spdobj1, spdobj2);
    reci(i) = spd_eqtest(spdobj1, spdobj2, "intrinsic");
end
        
%   print the results
sprintf(" empirical error with extrinsic dist : %.3f\n                      intrinsic dist : %.3f",(sum(recx<=0.05)+1)/(ntest+1),sum(reci<=0.05)/ntest)
