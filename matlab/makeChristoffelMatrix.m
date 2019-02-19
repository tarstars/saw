function [chm] = makeChristoffelMatrix(c, n)
chm = zeros(3, 3);
for p = 1:3
    for r = 1:3
        for q = 1:3
            for s = 1:3
                chm(p, r) = chm(p, r) + c(p, q, r, s) * n(q) * n(s);
            end
        end
    end
end