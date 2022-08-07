reset

set terminal epslatex standalone color size 5.,3.5 font ',17.28'
set output 'result.tex'

set xlabel '$t$'
set ylabel '$\omega_z$'

set xrange [0:10]
set yrange [0.1:1.1]

set xtics 0,5,10
set ytics 0.1,0.2,0.9

set format x '$% .0f$'
set format y '$% .1f$'

set style line 1 lc rgb '#FF0000'
set style line 2 lc rgb '#0000FF'
set style line 3 lc rgb '#33AA00'
set style line 4 lc rgb '#000000'

array filenames[4] = [ \
  '48/particle.dat', \
  '64/particle.dat', \
  '96/particle.dat', \
  'Aidun.dat', \
]

array titles[4] = [ \
  '$48$', \
  '$64$', \
  '$96$', \
  'ref.', \
]

set style line 5 lc rgb '#AAAAAA' lw 3
set key left top spacing 1.2 box ls 5

plot \
  filenames[1] u 1:7 t titles[1] ls 1 lw 5           w l, \
  filenames[2] u 1:7 t titles[2] ls 2 lw 5           w l, \
  filenames[3] u 1:7 t titles[3] ls 3 lw 5           w l, \
  filenames[4] u 1:2 t titles[4] ls 4 lw 3 pt 6 ps 2 w p

