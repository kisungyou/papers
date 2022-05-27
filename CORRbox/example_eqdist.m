%% EXAMPLE : test equality of distributions 
%            empirical Type 1 error at significance level 0.05.

%   preliminary
clear; close; clc; addpath(genpath(pwd));

%   example setup
ntest = 100;
p     = 5;              % size of correlation matrices
rec1  = zeros(1,ntest); % vector to record p-values via test2bg method
rec2  = zeros(1,ntest); %                               test2wass method

n1 = 10; % number of objects for class 1
n2 = 8;  %                       class 2

%   iterate !
for i=1:ntest
    dat1 = zeros(p,p,n1);
    dat2 = zeros(p,p,n2); 
    
    for m=1:n1
        dat1(:,:,m) = corr(randn(100,p));
    end
    for n=1:n2
        dat2(:,:,n) = corr(randn(100,p));
    end
    
    corr_obj1 = corr_initialize(dat1); % don't forget to initialize the data !
    corr_obj2 = corr_initialize(dat2); 
    
    rec1(i) = corr_test2bg(corr_obj1, corr_obj2);
    rec2(i) = corr_test2wass(corr_obj1, corr_obj2);
    fprintf("* iteration %d/%d complete\n",i,ntest);
end
        
%   print the results
sprintf("* empirical error with test2bg   : %.3f\n                       test2wass : %.3f",(sum(rec1<=0.05)+1)/(ntest+1),sum(rec2<=0.05)/ntest)
