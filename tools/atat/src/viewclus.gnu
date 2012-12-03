unset key
set term png
set out 'clus.png'
set size 0.75,1
set pointsize 10
unset xtics
unset ytics
unset ztics
splot 'clusl.tmp' w l 2,'clusp.tmp' w p 1 7
#pause -1
set term png
set out 'sh.png'
#set view ,,,1.5
set cbrange [-1.0:1.0]
unset colorbox
unset border
unset xtics
unset ytics
unset ztics
#set pm3d scansautomatic 
set pm3d at s
set pm3d scansbackward
set nosurf
splot 'sh.tmp' u 1:2:3:4 w l
#pause -1
