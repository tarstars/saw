function testMakeChristoffelMatrix()
c = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 6.6e10, 2.65e10);

range = 10.0 / 180 * pi;
angle_from = - range;
angle_to =  range;

points = 200;
phis = angle_from:((angle_to - angle_from)/points):angle_to;
x1=zeros(points);
y1=zeros(points);
x2=zeros(points);
y2=zeros(points);
x3=zeros(points);
y3=zeros(points);
ind = 0;
for phi = phis
    ind = ind + 1;
    n = [cos(phi), sin(phi), 0];
    chm = makeChristoffelMatrix(c, n);
    gammas = eig(chm);
    gammas = sort(gammas);
    t = 1;
    v = sqrt(gammas(t) / 5.96e3);
    x1(ind) = n(1) / v;
    y1(ind) = n(2) / v;
    t = 2;
    v = sqrt(gammas(t) / 5.96e3);
    x2(ind) = n(1) / v;
    y2(ind) = n(2) / v;
    t = 3;
    v = sqrt(gammas(t) / 5.96e3);
    x3(ind) = n(1) / v;
    y3(ind) = n(2) / v;
end
hold off
plot(x1, y1, '',  x2, y2, '', x3, y3, '')
hold on

points = 20;
phis = angle_from:((angle_to - angle_from)/points):angle_to;
for phi = phis
    n = [cos(phi), sin(phi), 0];
    chm = makeChristoffelMatrix(c, n);
    [q, gammas] = eig(chm.');
    for t = 1 : 3
        v = sqrt(gammas(t,t) / 5.96e3);
        x = n(1) / v;
        y = n(2) / v;
        px = q(1,t) * 0.25e-4;
        py = q(2,t) * 0.25e-4;
        quiver(x - px / 2, y - py / 2, px, py);
    end
end
