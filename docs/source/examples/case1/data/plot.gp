reset

set terminal epslatex standalone color size 5.,3.5 font ',17.28'
set output 'result.tex'

set xlabel '$t$'
set ylabel '$x$'

set xrange [0:150]
set yrange [0.25:0.55]

set xtics 0,50,150
set ytics 0.3,0.1,0.5

set format x '$% .0f$'
set format y '$% .1f$'

set style line 1 lc rgb '#FF0000'
set style line 2 lc rgb '#0000FF'
set style line 3 lc rgb '#33AA00'
set style line 4 lc rgb '#000000'

array filenames[4] = [ \
  '16/particle.dat', \
  '32/particle.dat', \
  '48/particle.dat', \
  'Feng1994.dat', \
]

array titles[4] = [ \
  '$16$', \
  '$32$', \
  '$48$', \
  'ref.', \
]

set style line 5 lc rgb '#AAAAAA' lw 3
set key right bottom spacing 1.2 box ls 5

plot \
  filenames[1] u 1:2 t titles[1] ls 1 lw 5           w l, \
  filenames[2] u 1:2 t titles[2] ls 2 lw 5           w l, \
  filenames[3] u 1:2 t titles[3] ls 3 lw 5           w l, \
  filenames[4] u 1:2 t titles[4] ls 4 lw 3 pt 6 ps 2 w p

