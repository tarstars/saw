#!/usr/bin/gnuplot

set terminal png
set output "mutant_z11.png"
set size square
set xrange [-0.002:0.002]
set yrange [-0.002:0.002]
plot "mutant_z11.txt"