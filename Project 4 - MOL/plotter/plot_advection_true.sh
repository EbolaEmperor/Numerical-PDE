echo "
set terminal postscript eps size 5,4
set xrange [-2:8]
set yrange [-0.4:1.1]
set output '../report/figures/advection-0.eps'
set sample 1000
plot exp(-20*((x-2)**2))+exp(-((x-5)**2)) w l title ''
exit" | gnuplot