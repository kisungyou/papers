function Sig = spdaux_2014CLX(p, snumber)

if (snumber==4)
    Omega = zeros(p,p);
    for i=1:(p-1)
        for j=(i+1):p
            theval = 0.8^abs(i-j);
            Omega(i,j) = theval;
            Omega(j,i) = theval;
        end
    end
    for i=1:p
        Omega(i,i) = 1;
    end
    Dhalf = diag(rand(p,1)*2 + 1);
    Sig   = Dhalf*(Omega\Dhalf);
elseif (snumber==3)
    Omega = zeros(p,p);
    for i=1:(p-1)
        Omega(i,i+1) = 0.8;
    end
    for i=1:(p-2)
        Omega(i,i+2) = 0.4;
    end
    for i=1:(p-3)
        Omega(i,i+3) = 0.2;
    end
    Omega = Omega + Omega';
    for i=1:p
        Omega(i,i) = 2;
    end
    Sig = (Omega\diag(ones(p,1)));
end
end
