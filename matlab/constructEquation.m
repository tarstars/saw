function [eq, mat] = constructEquation(c, rho)
s(1) = sym('s1');
s(2) = sym('s2');
s(3) = sym('s3');
gmat = sym('gmat', [3,3]);

for p = 1:3
    for q = 1:3
        gmat(p, q) = 0;        
        for r = 1:3
            for ss = 1:3
                gmat(p, q) = gmat(p, q) +  c(p, r, q, ss) * s(r) * s(ss); 
            end
        end
        if p == q
            gmat(p, q) = gmat(p, q) - rho;
        end
    end
end

mat = gmat;
eq = det(gmat);