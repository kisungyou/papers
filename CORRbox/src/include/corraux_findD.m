function D = corraux_findD(A,B)

% find D that minimizes d(A, DBD) with some default params
% see https://github.com/kisungyou/Riemann/blob/master/src/riemann_manifolds.h

C1 = A;
C2 = B;

p = size(A, 1);
maxiter = 50;
stopeps = 10^(-8);

Dold = eye(p);
Dnew = zeros(p,p);
Dtmp = zeros(p,p);
Ctmp = zeros(p,p);

Dhalf = zeros(p,p);
Dhalfinv = zeros(p,p);

Ddel1 = zeros(p,p);
Ddel2 = zeros(p,p);
Ddel3 = zeros(p,p);

costold = corraux_spddist(A,B);
costnew = 0;
Dincval = 0;

vec_delta = linspace(0.0001, 1.9999, 20);
vec_cost  = zeros(1,20);

for it=1:maxiter
    % compute intermediate values
    Ddel1 = Dold*real(logm(C2*Dold*(C1\Dold)));
    Ddel2 = 2.0*mat_symm(Ddel1, true);
    
    Dhalf    = mat_diaghalf(Dold);
    Dhalfinv = mat_diaginvhalf(Dold);
    Ddel3    = Dhalfinv*Ddel2*Dhalfinv;
    
    % iterate over values
    for i=1:20
        Dtmp = Dhalf*real(expm(-vec_delta(i)*Ddel3))*Dhalf;
        Ctmp = Dtmp*C2*Dtmp;
        vec_cost(i) = corraux_spddist(C1, Ctmp);
    end
    
    % optimal one
    [~,minid] = min(vec_cost);
    Dnew = Dhalf*real(expm(-vec_delta(minid)*Ddel3))*Dhalf;
    costnew = vec_cost(minid);
    
    % updating information
    Dincval = norm(Dold-Dnew,"fro")/norm(Dold,"fro");
    Dold = Dnew;
    
    cond1 = (Dincval < stopeps);
    cond2 = ((it > 1)&&(costnew>costold));
    if (cond1||cond2)
        break
    end
    costold = costnew;
end
D = Dold;
end


%% auxiliary functions
function B = mat_symm(A, use_diag)
    if (nargin < 2)
        use_diag = false;
    end
    B = (A+A')/2;
    if (use_diag==true)
        B = diag(diag(B));
    end
end

function output = mat_diaghalf(D)
    output = diag(sqrt(diag(D)));
end

function output = mat_diaginvhalf(D)
    output = diag(1./sqrt(diag(D)));
end
