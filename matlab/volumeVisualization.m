function volumeVisualization
n = 100;
a = zeros(n, n, n);
cx = n / 2;
cy = n / 2;
cz = n / 2;
sigma = 10;
gridx = a;
gridy = a;
gridz = a;
for p = 1:n
    for q = 1:n
        for r = 1:n
            gridx(p, q, r) = q;
            gridy(p, q, r) = p;
            gridz(p, q, r) = r;
            dx = p - cx;
            dy = q - cy;
            dz = r - cz;
            rad = 9 * dx * dx + dy * dy + dz * dz;
            rad = rad / (sigma * sigma);
            a(p, q, r) = exp( - rad);
        end
    end
end

fig1 = figure(1);
isosurface(gridx, gridy, gridz, a, 0.5);
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud

fig2 = figure(2);
b = fftn(a);
br = circshift(abs(b), [cx, cy, cz]);
isosurface(gridx, gridy, gridz, br, 0.5);
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
