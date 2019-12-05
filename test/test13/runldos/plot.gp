
set terminal qt size 800,1000 enhanced font "Helvetica,20" linewidth 2

ha2ev = 27.2114

set multiplot layout 2, 1

set xrange [-5:5]
set xlabel "Energy (eV)"
set ylabel "DOS" offset 0,0.5
p "ldos.out" u ($1*ha2ev):2 w l lw 1 lc rgb "red" title "LDOS"

unset xrange
set xlabel "z (bohr)"
set ylabel "Energy (eV)"
p "band_edges.out" u 1:($2*ha2ev) w lp lw 2 pt 7 lc rgb "red" title "VBM", \
   "" u 1:($3*ha2ev) w lp lw 2 pt 7 lc rgb "blue" title "CBm"

unset multiplot
unset terminal