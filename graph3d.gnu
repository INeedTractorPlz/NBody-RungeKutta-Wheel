#! /usr/bin/gnuplot -persist


set xrange [-2:2]
set yrange [-2:2]
set zrange [-0.0025:0.003]

splot "NBodyLP_0.dat" w l,"NBodyLP_1.dat" w l,"NBodyLP_2.dat" w l
pause -1
