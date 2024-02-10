reset

load 'fit-circle.dat'

min(x, y) = x > y ? y : x
max(x, y) = x < y ? y : x

xe(a, t) = a*cos(t)
ye(b, t) = b*sin(t)

x_t0 = atan2(-e0_b*sin(e0_theta), e0_a*cos(e0_theta))
x_t1 = x_t0+pi
y_t0 = atan2(+e0_b*cos(e0_theta), e0_a*sin(e0_theta))
y_t1 = y_t0+pi
x0 = e0_x+xe(e0_a, x_t0)*cos(e0_theta)-ye(e0_b, x_t0)*sin(e0_theta)
x1 = e0_x+xe(e0_a, x_t1)*cos(e0_theta)-ye(e0_b, x_t1)*sin(e0_theta)
y0 = e0_y+xe(e0_a, y_t0)*sin(e0_theta)+ye(e0_b, y_t0)*cos(e0_theta)
y1 = e0_y+xe(e0_a, y_t1)*sin(e0_theta)+ye(e0_b, y_t1)*cos(e0_theta)
xmin = min(e0_xp, min(x0, x1))
xmax = max(e0_xp, max(x0, x1))
ymin = min(e0_yp, min(y0, y1))
ymax = max(e0_yp, max(y0, y1))

mrgn = 0.25

xmin = xmin-mrgn
xmax = xmax+mrgn
ymin = ymin-mrgn
ymax = ymax+mrgn

lx = xmax-xmin
ly = ymax-ymin

set terminal epslatex standalone color size lx,ly font ',17.28'
set output 'fit-circle.tex'

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
e0x = e0_x
e0y = e0_y
e1x = e0_x+e0_a*cos(e0_theta)
e1y = e0_y+e0_a*sin(e0_theta)
e2x = e0_x-e0_b*sin(e0_theta)
e2y = e0_y+e0_b*cos(e0_theta)

set object circle at first e0_xc, first e0_yc size first e0_r fillstyle empty fc rgb '#0000FF' lw 3

set arrow from first e0_xp, first e0_yp to first e0_xc, first e0_yc as 1
set label '$\left( x_c, y_c \right)$' right at first e0_xc, first e0_yc
set label '$\left( x_p, y_p \right)$\ \ ' right at first e0_xp, first e0_yp

plot \
  NaN notitle
