set terminal png
set output "iter.png"
set xlabel "h"
set ylabel "iteration count"
plot './data/iter.dat' using 1:2 title ""