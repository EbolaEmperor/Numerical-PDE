./../heat ../tests/$1.json
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [$2:$3]
set output '../report/figures/$1-1.eps'
plot 'result1.txt' w lp pt 7 title ''
exit
" | gnuplot
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [$4:$5]
set output '../report/figures/$1-2.eps'
plot 'result2.txt' w lp pt 7 title ''
exit
" | gnuplot
echo "
set terminal postscript eps size 5,4
set xrange [0:1]
set yrange [$6:$7]
set output '../report/figures/$1-3.eps'
plot 'result10.txt' w lp pt 7 title ''
exit
" | gnuplot