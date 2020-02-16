#!/usr/bin/gnuplot
reset
set term pdf enhanced font 'Helvetica,18'
set output 'matrix.pdf'
set key outside
set tics in scale 0.5
set size ratio -1
set ytics 10
set xtics 10
set mxtics 2
set mytics 2
set tics out nomirror
set datafile separator " "

set style fill solid 1.0 border rgb 'black'
set style line 2604 linetype -1 linewidth .7
set colorbox border 2604
set cbtics .2

set palette maxcolors 10000
set palette defined (0  "red", 0.5 "grey", 1 "blue")

set xlabel 'Cα';
set ylabel 'Cα';

plot 'text-output.xvg' using 1:2:(2.8):3 notitle w points pt 7 ps .5 lc palette;
