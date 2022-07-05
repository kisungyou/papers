function dval = corraux_wassersteinD(D, p, wx, wy)

% https://github.com/kisungyou/Riemann -> auxiliary_ported
if (nargin < 2)
    p = 2.0;
end
dxy   = D;
p     = max(1, double(p));
[m,n] = size(D);

if (nargin < 3)
    wx = ones(m,1)/m;
end
if (nargin < 4)
    wy = ones(n,1)/n;
end

wx  = wx/sum(wx); wx = wx(:); % vector is a single-column matrix
wy  = wy/sum(wy); wy = wy(:); 

cxy   = dxy.^p;

c  = cxy(:);
A1 = kron(ones(1,n), eye(m));
A2 = kron(eye(n), ones(1,m));
A  = [A1;A2];

f_obj = c;
f_con = A;
f_rhs = [wx;wy]; f_rhs = f_rhs(:);
f_lwb = zeros(length(f_obj),1);
f_upb = ones(length(f_obj),1);
f_opt = optimoptions("linprog","Display","none");
f_sol = linprog(f_obj, [], [], f_con, f_rhs, f_lwb, f_upb, f_opt);

dval  = (sum(f_sol.*c))^(1/p);

end