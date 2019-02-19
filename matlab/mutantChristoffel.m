function [s3p, q] = mutantChristoffel(mutantPoly, mutantMatrix, sx, sy)

mutantPoly1 = subs(mutantPoly, 's1', 0);
mutantMatrix1 = subs(mutantMatrix, 's1', 0);


mutantPoly2 = subs(mutantPoly1, 's2', sy);
mutantMatrix2 = subs(mutantMatrix1, 's2', sy);

s3p = double(solve(mutantPoly2, 's3'));
s3p = s3p.';

for s3c = s3p
    valueOfPoly = double(subs(mutantPoly2, 's3', (s3c)))
    zdetMatrix = subs(mutantMatrix2, 's3', (s3c));
    determ = det(zdetMatrix)
end
q = []