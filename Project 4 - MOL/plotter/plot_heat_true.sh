./../trueSol 0 1e-3
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [-0.001:1]
set output '../report/figures/heat-a-0-0.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./../trueSol 0.0025 1e-3
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [-0.001:0.5]
set output '../report/figures/heat-a-0-1.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./../trueSol 0.005 1e-3
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [-0.001:0.5]
set output '../report/figures/heat-a-0-2.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./../trueSol 0.025 1e-3
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [-0.001:0.5]
set output '../report/figures/heat-a-0-3.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot