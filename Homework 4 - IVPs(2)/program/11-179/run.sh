#!bin/bash

./test
echo "
set terminal postscript eps color
set output '11-179-1.eps'
plot 'Backward Eular' w lp dashtype 2 pt 7, 0.5*exp(-1e6*x)+cos(x)
exit
" | gnuplot
echo "
set terminal postscript eps color
set output '11-179-2.eps'
plot 'Trapezoidal' w lp dashtype 2 pt 7, 0.5*exp(-1e6*x)+cos(x)
exit
" | gnuplot