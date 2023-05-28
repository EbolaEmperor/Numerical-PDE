#!/bin/bash

./solve samples/orbit1.json
echo "
set terminal eps size 5,4
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/0-1.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve samples/orbit2.json
echo "
set terminal eps size 5,4
set output 'report/figures/0-2.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve samples/orbit3.json
echo "
set terminal eps size 5,4
set xrange [-1:1]
set output 'report/figures/0-3.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve samples/orbit4.json
echo "
set terminal eps size 5,4
set output 'report/figures/0-4.eps'
set ticslevel 0
set view 60,30,1,1.5
splot 'result-dense-discrete.txt' u 2:3:4 w l title ''
exit
" | gnuplot
./solve samples/orbit5.json
echo "
set terminal eps size 8,6
set output 'report/figures/0-5.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title '', 'result-dense-discrete.txt' u 5:6 w l title '', 'result-dense-discrete.txt' u 8:9 w l title ''
exit
" | gnuplot
./solve samples/orbit6.json
echo "
set terminal eps size 8,6
set output 'report/figures/0-6.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title '', 'result-dense-discrete.txt' u 5:6 w l title '', 'result-dense-discrete.txt' u 8:9 w l title ''
exit
" | gnuplot
./solve tests/0-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/0-7.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./solve tests/1-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/1-1.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/1-2-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/1-2.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/2-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/2-1.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/2-2-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/2-2.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/3-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/3-1.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/3-3-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/3-3.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/3-5-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/3-5.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/3-2-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/3-2.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/3-4-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/3-4.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/3-6-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/3-6.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/3-7-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/3-7.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./solve tests/3-7-2.json
echo "
set terminal eps size 8,6
set output 'report/figures/3-8.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./solve tests/3-7-3.json
echo "
set terminal eps size 8,6
set output 'report/figures/3-9.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./solve tests/4-1-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/4-1.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/4-2-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/4-2.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/4-3-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/4-3.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/4-4-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/4-4.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/5-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/5-1.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/5-2-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/5-2.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/5-3-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/5-3.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/5-4-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/5-4.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/5-6-1.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/5-6.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/5-6-2.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/5-7.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/5-5-1.json
echo "
set terminal png size 800,600
set output 'report/figures/5-5.png'
plot 'result.txt' w lp pt 6 ps 1 title ''
exit
" | gnuplot
./solve tests/6-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/6-1.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/6-2-5.json
echo "
set terminal eps size 8,6
set output 'report/figures/6-2.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/7-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/7-1.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/7-2-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/7-2.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/7-3-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/7-3.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/7-4-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/7-4.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/7-5-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/7-5.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/7-6-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/7-6.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/7-7-1.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/7-7.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/7-7-2.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/7-8.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/7-7-3.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/7-9.eps'
plot 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/8-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/8-1.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/8-2-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/8-2.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/9-1-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/9-1.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/9-2-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/9-2.eps'
plot 'result.txt' u 2:3 w p pt 6 title '', 'result-dense-discrete.txt' u 2:3 w l lt 1 title ''
exit
" | gnuplot
./solve tests/9-3-1.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/9-3.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/9-3-2.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/9-4.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/9-3-3.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/9-5.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/9-4-1.json
echo "
set terminal eps size 8,6
set yrange [-2.5:2.5]
set output 'report/figures/9-6.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/9-5-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/9-7.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/9-6-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/9-8.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/10-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/10-1.eps'
plot 'result.txt' w l title ''
exit
" | gnuplot
./solve tests/10-2-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/10-2.eps'
plot 'result.txt' w lp pt 6 title ''
exit
" | gnuplot
./solve tests/10-3-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/10-3.eps'
plot 'result.txt' w p pt 6 title '', 'result-dense-discrete.txt' w l lt 1 title ''
exit
" | gnuplot
./solve tests/10-4-1.json
echo "
set terminal eps size 8,6
set output 'report/figures/10-4.eps'
plot 'result.txt' w p pt 6 title '', 'result-dense-discrete.txt' w l lt 1 title ''
exit
" | gnuplot
./solve tests/11-1.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/11-2.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-2-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/11-3.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-3-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/11-4.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-4-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/11-5.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-5-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/11-6.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-6-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/11-7.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-7-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/11-8.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-8-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/11-9.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/11-9-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/11-10.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-1-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/12-1.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-2-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/12-2.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-3-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/12-3.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-4-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/12-4.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-5-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/12-5.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-6-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/12-6.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-7-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/12-7.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/12-8-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/12-8.eps'
plot 'result.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-1-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/13-1.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-2-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/13-2.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-3-5.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/13-3.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-4-4.json
echo "
set terminal eps size 8,6
set xrange [-1.4:1.1]
set yrange [-1.3:1.3]
set output 'report/figures/13-4.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-5-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/13-5.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-6-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/13-6.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-7-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/13-7.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/13-8-4.json
echo "
set terminal eps size 8,6
set output 'report/figures/13-8.eps'
plot 'result-dense-discrete.txt' u 2:3 w l title ''
exit
" | gnuplot
./solve tests/14-1-1.json
mv result.txt result1.txt
./solve tests/14-1-2.json
mv result.txt result2.txt
./solve tests/14-1-3.json
mv result-dense-discrete.txt result3.txt
echo "
set terminal eps size 8,6
set output 'report/figures/14-1.eps'
plot 'result1.txt' u 2:3 w l title 'Eular, 24000 Steps', 'result2.txt' u 2:3 w l title 'Classical RK, 6000 Steps', 'result3.txt' u 2:3 w l lt 4 title 'DOPRI5, 100 Steps'
exit
" | gnuplot
