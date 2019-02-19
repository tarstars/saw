data = importdata('te_biss_z.dat');
axis equal
axis square
x=data(:, 1);
y=data(:, 2);
scatter(x,y,2)