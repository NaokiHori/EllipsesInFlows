np = system('find ./output/log -type f -name particle*.dat | wc -l')

set yrange [0.:1.]

plot \
  for[i=0:np-1:1] \
  sprintf('output/log/particle%010d.dat', i) u 1:2 notitle w l
