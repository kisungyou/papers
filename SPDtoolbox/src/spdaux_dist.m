function dval = spdaux_dist(X,Y,mode)

% distance between two SPD matrices
if (nargin<3)
    mode = 'extrinsic';
end

if strcmp(mode,'intrinsic')
    trnorm = @(A) sqrt(trace(A*A));
    dval   = real(trnorm(real(logm(X\Y)))); % AIRM/LERM; not clear what to use.
else
    p  = size(X,1);
    xx = reshape(real(logm(X)),[p^2,1]);
    yy = reshape(real(logm(Y)),[p^2,1]);
    dval = sqrt(sum((xx-yy).^2));
end
end