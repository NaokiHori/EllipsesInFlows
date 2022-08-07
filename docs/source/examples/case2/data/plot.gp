reset

set terminal epslatex standalone color size 5.,3.5 font ',17.28'
set output 'result.tex'

set xlabel '$y$'
set ylabel '$x$'

set xrange [0:25]
set yrange [0.26:0.42]

set xtics 0,5,25
set ytics 0.26,0.04,0.42

set format x '$% .0f$'
set format y '$% .2f$'

set style line 1 lc rgb '#FF0000'
set style line 2 lc rgb '#0000FF'
set style line 3 lc rgb '#33AA00'
set style line 4 lc rgb '#000000'

array filenames[4] = [ \
  '32/particle.dat', \
  '64/particle.dat', \
  '96/particle.dat', \
  'PanGlowinski2002.dat', \
]

array titles[4] = [ \
  '$32$', \
  '$64$', \
  '$96$', \
  'ref.', \
]

set style line 5 lc rgb '#AAAAAA' lw 3
set key right top spacing 1.2 box ls 5

plot \
  filenames[1] u 2:1 t titles[1] ls 1 lw 5           w l, \
  filenames[2] u 2:1 t titles[2] ls 2 lw 5           w l, \
  filenames[3] u 2:1 t titles[3] ls 3 lw 5           w l, \
  filenames[4] u 2:1 t titles[4] ls 4 lw 3 pt 6 ps 2 w p

