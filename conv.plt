set terminal png
set output "conv.png"
set xlabel "n"
set ylabel "maximum norm of error"
plot './data/conv.dat' using 1:2 title ""
