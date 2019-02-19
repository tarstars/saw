data = importdata('te_volume.dat');
axis equal
axis square
x=data(:, 1);
y=data(:, 2);
z=data(:, 3);
scatter3(x,y,z)