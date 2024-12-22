set term x11 persist size 1200,800
set multiplot layout 3,2 title 'Wave Equation Solutions'
set grid xtics ytics ls 3
set grid lt 1 lc rgb '#dddddd'
set style line 1 lc rgb '#0000FF' lt 1 lw 2
set style line 2 lc rgb '#FF0000' lt 1 lw 2
set style line 3 lc rgb '#00FF00' lt 1 lw 2
set xtics 0.2
set ytics 0.2
set xrange [0:1]
set title 'First-order Numerical Solution'
set ylabel 't'
set yrange [0:1]
set view map
set palette defined (0 '#000080', 0.5 '#00FF00', 1 '#FFFF00')
set cbrange [0:1.6]
plot 'first_order_numerical.dat' using 1:2:3 with image notitle
set title 'Second-order Numerical Solution'
plot 'second_order_numerical.dat' using 1:2:3 with image notitle
set title 'First-order Error'
set palette defined (0 '#000080', 0.5 '#00FF00', 1 '#FFFF00')
set cbrange [0:0.35]
plot 'first_order_error.dat' using 1:2:3 with image notitle
set title 'Second-order Error'
plot 'second_order_error.dat' using 1:2:3 with image notitle
unset view
set origin 0.0,0.0
set size 1.0,0.34
set title 'Second-order solution at t = 0.50'
set ylabel 'u(x,t)'
set yrange [0.4:1.3]
set key bottom right
plot 'first_order_comparison.dat' using 1:2 title 'First-order' with lines ls 1, \
     'second_order_comparison.dat' using 1:2 title 'Second-order' with lines ls 2, \
     'first_order_comparison.dat' using 1:3 title 'Exact' with lines dashtype 2 lc rgb '#FF0000'
unset multiplot
