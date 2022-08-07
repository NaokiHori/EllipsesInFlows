beta = 2.

f0(x) = 0.5     *(1.+tanh(x)    )
f1(x) = 0.5*beta*(1.-tanh(x)**2.)

reset
set terminal epslatex standalone color size 5.,3.5 font ',12'
set output 'indicator.tex'
set xlabel '$x$'
set ylabel '$f \left( x \right), f^{\prime} \left( x \right)$'
xmin = -4.
xmax = +4.
set xrange [xmin:xmax]
set yrange [-0.1:+1.1]
# set [xy]tics <start>, <incr>, <end>
set xtics xmin, 1., xmax
set ytics 0.5
set format x '$% .0f$'
set format y '$% .1f$'
set style line 1 lc rgb '#FF0000' lw 5
set style line 2 lc rgb '#0000FF' lw 5
set key at graph 0.3, graph 0.9 spacing 1.2
set samples 250
plot \
  f0(x) title '$f          \left( x \right)$' ls 1, \
  f1(x) title '$f^{\prime} \left( x \right)$' ls 2

