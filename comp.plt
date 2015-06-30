set terminal png
set output "comp.png"
set xlabel "n"
set ylabel "maximum norm of residual"
plot './data/comp.dat' using 1:2 title ""
