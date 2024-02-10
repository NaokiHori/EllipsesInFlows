
do for [case = 0:1:1] {

  reset

  load sprintf('fit-circles-%d.dat', case)

  array e_xs[2] = [e0_x, e1_x]
  array e_ys[2] = [e0_y, e1_y]
  array e_as[2] = [e0_a, e1_a]
  array e_bs[2] = [e0_b, e1_b]
  array e_ts[2] = [e0_theta, e1_theta]
  array e_xcs[2] = [e0_xc, e1_xc]
  array e_ycs[2] = [e0_yc, e1_yc]
  array e_rs[2] = [e0_r, e1_r]

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
    x2 = e_xcs[n]-e_rs[n]
    xmin = min(xmin, min(e_xcs[n]-e_rs[n], min(x0, x1)))
    xmax = max(xmax, max(e_xcs[n]+e_rs[n], max(x0, x1)))
    ymin = min(ymin, min(e_ycs[n]-e_rs[n], min(y0, y1)))
    ymax = max(ymax, max(e_ycs[n]+e_rs[n], max(y0, y1)))
  }

  mrgn = 0.25

  xmin = xmin-mrgn
  xmax = xmax+mrgn
  ymin = ymin-mrgn
  ymax = ymax+mrgn

  lx = xmax-xmin
  ly = ymax-ymin

  set terminal epslatex standalone color size lx,ly font ',17.28'
  set output sprintf('fit-circles-%d.tex', case)

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

  set style line 1 lc rgb '#000000' lw 7
  set style arrow 1 heads  ls 1
  set style arrow 2 nohead ls 1 dt 3
  asize = 0.5

  set object ellipse at first e0_x, first e0_y size first 2.*e0_a, first 2.*e0_b angle 180./pi*e0_theta fillstyle empty fc rgb '#FF0000' lw 3
  set object circle at first e0_xc, first e0_yc size first e0_r fillstyle empty fc rgb '#0000FF' lw 3
  set object ellipse at first e1_x, first e1_y size first 2.*e1_a, first 2.*e1_b angle 180./pi*e1_theta fillstyle empty fc rgb '#FF0000' lw 3
  set object circle at first e1_xc, first e1_yc size first e1_r fillstyle empty fc rgb '#0000FF' lw 3

  plot \
    NaN notitle
}

