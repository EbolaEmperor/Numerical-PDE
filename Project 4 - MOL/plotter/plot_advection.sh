./../advection ../tests/$1.json
echo "
set terminal postscript eps size 5,4
set xrange [15:25]
set yrange [-0.4:1.1]
set output '../report/figures/$1.eps'
set sample 1000
plot 'result.txt' w lp pt 7 ps 0.4 title '', exp(-20*((x-2-17)**2))+exp(-((x-5-17)**2)) w l dashtype '.' title ''
exit" | gnuplot