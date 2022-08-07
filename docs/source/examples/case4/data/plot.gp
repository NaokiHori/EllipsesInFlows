reset

set terminal epslatex standalone color size 5.,3.5 font ',17.28'
set output 'vfrac.tex'

set xlabel '$x$'
set ylabel '$\left\langle \phi \right\rangle_{y,t}$'

set xrange [0.0:1.0]
set yrange [0.0:0.5]

set xtics 0.00,0.25,1.00
set ytics 0.00,0.10,0.50

set style line 1 lc rgb '#FF0000'
set style line 2 lc rgb '#0000FF'
set style line 3 lc rgb '#000000'
set style line 4 lc rgb '#AAAAAA'

set key left top spacing 1.2 box ls 4 lw 2

plot \
  'vfrac.dat'    u 1:2 t 'result' ls 1 lw 2 w l, \
  'Chen2012.dat' u 2:1 t 'ref.  ' ls 2 pt 6 ps 1 w p, \
  0.167 notitle ls 3 lw 2 dt 2 w l

