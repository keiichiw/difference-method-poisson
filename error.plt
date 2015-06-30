set terminal png
set output "/dev/null"
set xlabel "log(h)"
set ylabel "log (error norm)"
set logscale x
set logscale y
plot './data/error.dat' using 1:2 title ""
set output "error.png"
replot 1/x title "y = 1/x"