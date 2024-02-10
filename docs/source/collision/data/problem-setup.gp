reset

array e_xs[2] = [0.5, 2.5]
array e_ys[2] = [0.5, 3.0]
array e_as[2] = [2.0, 1.5]
array e_bs[2] = [1.5, 1.0]
array e_ts[2] = [0.5, 1.2]

min(x, y) = x > y ? y : x
max(x, y) = x < y ? y : x

xe(a, t) = a*cos(t)
ye(b, t) = b*sin(t)

xmin = +1.e+16
xmax = -1.e+16
ymin = +1.e+16
ymax = -1.e+16

do for [n = 1:2:1] {
  x_t0 = atan2(-e_bs[n]*sin(e_ts[n]), e_as[n]*cos(e_ts[n]))
  x_t1 = x_t0+pi
  y_t0 = atan2(+e_bs[n]*cos(e_ts[n]), e_as[n]*sin(e_ts[n]))
  y_t1 = y_t0+pi
  x0 = e_xs[n]+xe(e_as[n], x_t0)*cos(e_ts[n])-ye(e_bs[n], x_t0)*sin(e_ts[n])
  x1 = e_xs[n]+xe(e_as[n], x_t1)*cos(e_ts[n])-ye(e_bs[n], x_t1)*sin(e_ts[n])
  y0 = e_ys[n]+xe(e_as[n], y_t0)*sin(e_ts[n])+ye(e_bs[n], y_t0)*cos(e_ts[n])
  y1 = e_ys[n]+xe(e_as[n], y_t1)*sin(e_ts[n])+ye(e_bs[n], y_t1)*cos(e_ts[n])
  xmin = min(xmin, min(x0, x1))
  xmax = max(xmax, max(x0, x1))
  ymin = min(ymin, min(y0, y1))
  ymax = max(ymax, max(y0, y1))
}

mrgn = 0.25

xmin = xmin-mrgn
xmax = xmax+mrgn
ymin = ymin-mrgn
ymax = ymax+mrgn

lx = xmax-xmin
ly = ymax-ymin

set terminal epslatex standalone color size lx,ly font ',17.28'
set output 'problem-setup.tex'

unset border

set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

unset xlabel
unset ylabel

set xrange [xmin:xmax]
set yrange [ymin:ymax]

unset xtics
unset ytics

set size ratio -1

set style line 1 lc rgb '#AAAAAA' lw 5
set style arrow 1 head   ls 1
set style arrow 2 nohead ls 1 dt 3
asize = 0.5

do for [n = 1 : 2 : 1] {
  set object ellipse at first e_xs[n], first e_ys[n] size first 2.*e_as[n], first 2.*e_bs[n] angle 180./pi*e_ts[n] fillstyle empty fc rgb '#FF0000' lw 3
  e0x = e_xs[n]
  e0y = e_ys[n]
  e1x = e_xs[n]+e_as[n]*cos(e_ts[n])
  e1y = e_ys[n]+e_as[n]*sin(e_ts[n])
  e2x = e_xs[n]-e_bs[n]*sin(e_ts[n])
  e2y = e_ys[n]+e_bs[n]*cos(e_ts[n])
  set arrow from first e0x, first e0y to first e1x, first e1y as 1
  set arrow from first e0x, first e0y to first e2x, first e2y as 1
  set arrow from graph 0.,  first e0y to graph 1.,  first e0y as 2
  set object circle at first e0x, first e0y size first 0.5 arc [0.:180./pi*e_ts[n]] fillstyle empty fc rgb '#AAAAAA' lw 3
  set label sprintf('$a_%d$', n-1) center at first 0.5*(e0x+e1x), first 0.5*(e0y+e1y) front
  set label sprintf('$b_%d$', n-1) center at first 0.5*(e0x+e2x), first 0.5*(e0y+e2y) front
  set label sprintf('$\theta_%d$', n-1) center at first e0x+0.5*cos(0.5*e_ts[n]), first e0y+0.5*sin(0.5*e_ts[n]) front
  set label sprintf('$( x_%d, y_%d )$', n-1, n-1) center at first e0x, first e0y-0.25 front
}

plot \
  NaN notitle
