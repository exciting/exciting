

! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program main

use inputdom
use modinput
implicit none
! read the input parameters from spacegroup.in

call loadinputDOM("spacegroup.xml")
input=getstructinput(inputnp)
call  initatomcounters()
! generate the lattice and atomic basis
call gencrystal
! write the structural data to GEOMETRY.OUT

! write the XCrySDen and V_Sim files for plotting
call writegeomxml
stop
end program

!BOI
! !TITLE: The Spacegroup Manual\\ Version 1.1.2
! !AUTHORS: J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
! {\sf Spacegroup} is a utility which produces crystal geometry for use with the
! {\sf EXCITING} code, from the space group defined by its Hermann-Mauguin
! symbol and lattice vector lengths and angles. {\sf Spacegroup} recognises all
! 230 space groups in various coordinate settings giving a total of 530 possible
! symbols, which are tabulated below. The code also provides output compatible
! with the {\sf XCrysDen} or {\sf V\_Sim} packages for visualisation of the
! crystal structure.
! \section{Usage}
! Only one input file, {\tt spacegroup.in}, is required. The structure of this
! file is illustrated by the following example for the high $T_c$ superconductor
! La$_2$CuO$_4$:
! \begin{verbatim}
!  'Bmab'                                : hrmg
!  10.0605232 10.0605232  24.972729      : a, b, c
!   90.0     90.0     90.0               : ab, ac, bc
!   1  1  1                              : ncell
!   .true.                               : primcell
!   3                                    : nspecies
!   'La' 'La.in'                         : spsymb, spfname
!   1                                    : nwpos
!   0.0000    0.0000    0.3608           : wpos
!   'Cu' 'Cu.in'
!   1
!   0.0000    0.0000    0.0000
!   'O' 'O.in'
!   2
!   0.2500    0.2500    0.0000
!   0.0000    0.0000    0.1820
! \end{verbatim}
! The input parameters are defined as follows:
! \vskip 6pt
! {\tt hrmg}\\
! The Hermann-Mauguin symbol of a space group listed in the table below
! (case-sensitive)
! \vskip 6pt
! {\tt a}, {\tt b}, {\tt c}\\
! Lattice vector lengths in Bohr (i.e. atomic units, {\bf NOT} \AA ngstroms)
! \vskip 6pt
! {\tt ab}, {\tt ac}, {\tt bc}\\
! Angles in degrees between lattice vectors ${\bf a}$ and ${\bf b}$; ${\bf a}$
! and ${\bf c}$; and ${\bf b}$ and ${\bf c}$, respectively
! \vskip 6pt
! {\tt ncell}\\
! The number of unit cells required in each direction
! \vskip 6pt
! {\tt primcell}\\
! Set to {\tt .true.} if the primitive unit cell should be found
! \vskip 6pt
! {\tt nspecies}\\
! Number of atomic species
! \vskip 6pt
! {\tt spsymb}, {\tt spfname}\\
! The species symbol and the species filename required by {\sf EXCITING}
! \vskip 6pt
! {\tt nwpos}\\
! The number of Wyckoff positional coordinates
! \vskip 6pt
! {\tt wpos}\\
! Wyckoff positional coordinates in fractions of the lattice vectors
! \vskip 6pt
! Note that {\tt nwpos} and {\tt wpos} are repeated as many times as there are
! species. After creating the input file, the {\tt spacegroup} command is run
! and the files {\tt GEOMETRY.OUT} and {\tt crystal.xsf} should be produced.
! The {\tt GEOMETRY.OUT} file can simply be appended to an {\tt exciting.in}
! file. If {\sf XCrysDen} is available, then use the command
! \vskip 6pt
! \hskip 24pt {\tt xcrysden --xsf crystal.xsf}
! \vskip 6pt
! to render the unit cell.
! \section{Table of space group symbols}
! We acknowledge Ralf W. Grosse-Kunstleve ({\tt http://cci.lbl.gov/sginfo/})
! for the following table which associates space group numbers, Sch\"{o}nflies
! symbols, Hermann-Mauguin symbols, and Hall symbols.
! \newpage
! \begin{center}
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!    1        &  C$_{1}^{1}$   &  P1                          &   P 1      \\
!    2        &  C$_{i}^{1}$   &  P-1                         &  -P 1      \\
!    3:b      &  C$_{2}^{1}$   &  P2:b = P121                 &   P 2y     \\
!    3:c      &  C$_{2}^{1}$   &  P2:c = P112                 &   P 2      \\
!    3:a      &  C$_{2}^{1}$   &  P2:a = P211                 &   P 2x     \\
!    4:b      &  C$_{2}^{2}$   &  P21:b = P1211               &   P 2yb    \\
!    4:c      &  C$_{2}^{2}$   &  P21:c = P1121               &   P 2c     \\
!    4:a      &  C$_{2}^{2}$   &  P21:a = P2111               &   P 2xa    \\
!    5:b1     &  C$_{2}^{3}$   &  C2:b1 = C121                &   C 2y     \\
!    5:b2     &  C$_{2}^{3}$   &  C2:b2 = A121                &   A 2y     \\
!    5:b3     &  C$_{2}^{3}$   &  C2:b3 = I121                &   I 2y     \\
!    5:c1     &  C$_{2}^{3}$   &  C2:c1 = A112                &   A 2      \\
!    5:c2     &  C$_{2}^{3}$   &  C2:c2 = B112 = B2           &   B 2      \\
!    5:c3     &  C$_{2}^{3}$   &  C2:c3 = I112                &   I 2      \\
!    5:a1     &  C$_{2}^{3}$   &  C2:a1 = B211                &   B 2x     \\
!    5:a2     &  C$_{2}^{3}$   &  C2:a2 = C211                &   C 2x     \\
!    5:a3     &  C$_{2}^{3}$   &  C2:a3 = I211                &   I 2x     \\
!    6:b      &  C$_{s}^{1}$   &  Pm:b = P1m1                 &   P -2y    \\
!    6:c      &  C$_{s}^{1}$   &  Pm:c = P11m                 &   P -2     \\
!    6:a      &  C$_{s}^{1}$   &  Pm:a = Pm11                 &   P -2x    \\
!    7:b1     &  C$_{s}^{2}$   &  Pc:b1 = P1c1                &   P -2yc   \\
!    7:b2     &  C$_{s}^{2}$   &  Pc:b2 = P1n1                &   P -2yac  \\
!    7:b3     &  C$_{s}^{2}$   &  Pc:b3 = P1a1                &   P -2ya   \\
!    7:c1     &  C$_{s}^{2}$   &  Pc:c1 = P11a                &   P -2a    \\
!    7:c2     &  C$_{s}^{2}$   &  Pc:c2 = P11n                &   P -2ab   \\
!    7:c3     &  C$_{s}^{2}$   &  Pc:c3 = P11b = Pb           &   P -2b    \\
!    7:a1     &  C$_{s}^{2}$   &  Pc:a1 = Pb11                &   P -2xb   \\
!    7:a2     &  C$_{s}^{2}$   &  Pc:a2 = Pn11                &   P -2xbc  \\
!    7:a3     &  C$_{s}^{2}$   &  Pc:a3 = Pc11                &   P -2xc   \\
!    8:b1     &  C$_{s}^{3}$   &  Cm:b1 = C1m1                &   C -2y    \\
!    8:b2     &  C$_{s}^{3}$   &  Cm:b2 = A1m1                &   A -2y    \\
!    8:b3     &  C$_{s}^{3}$   &  Cm:b3 = I1m1                &   I -2y    \\
!    8:c1     &  C$_{s}^{3}$   &  Cm:c1 = A11m                &   A -2     \\
!    8:c2     &  C$_{s}^{3}$   &  Cm:c2 = B11m = Bm           &   B -2     \\
!    8:c3     &  C$_{s}^{3}$   &  Cm:c3 = I11m                &   I -2     \\
!    8:a1     &  C$_{s}^{3}$   &  Cm:a1 = Bm11                &   B -2x    \\
!    8:a2     &  C$_{s}^{3}$   &  Cm:a2 = Cm11                &   C -2x    \\
!    8:a3     &  C$_{s}^{3}$   &  Cm:a3 = Im11                &   I -2x    \\
!    9:b1     &  C$_{s}^{4}$   &  Cc:b1 = C1c1                &   C -2yc   \\
!    9:b2     &  C$_{s}^{4}$   &  Cc:b2 = A1n1                &   A -2yac  \\
!    9:b3     &  C$_{s}^{4}$   &  Cc:b3 = I1a1                &   I -2ya   \\
!    9:-b1    &  C$_{s}^{4}$   &  Cc:-b1 = A1a1               &   A -2ya   \\
!    9:-b2    &  C$_{s}^{4}$   &  Cc:-b2 = C1n1               &   C -2ybc  \\
!    9:-b3    &  C$_{s}^{4}$   &  Cc:-b3 = I1c1               &   I -2yc   \\
!    9:c1     &  C$_{s}^{4}$   &  Cc:c1 = A11a                &   A -2a    \\
!    9:c2     &  C$_{s}^{4}$   &  Cc:c2 = B11n                &   B -2bc   \\
!    9:c3     &  C$_{s}^{4}$   &  Cc:c3 = I11b                &   I -2b    \\
!    9:-c1    &  C$_{s}^{4}$   &  Cc:-c1 = B11b = Bb          &   B -2b    \\
!    9:-c2    &  C$_{s}^{4}$   &  Cc:-c2 = A11n               &   A -2ac   \\
!    9:-c3    &  C$_{s}^{4}$   &  Cc:-c3 = I11a               &   I -2a    \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!    9:a1     &  C$_{s}^{4}$   &  Cc:a1 = Bb11                &   B -2xb   \\
!    9:a2     &  C$_{s}^{4}$   &  Cc:a2 = Cn11                &   C -2xbc  \\
!    9:a3     &  C$_{s}^{4}$   &  Cc:a3 = Ic11                &   I -2xc   \\
!    9:-a1    &  C$_{s}^{4}$   &  Cc:-a1 = Cc11               &   C -2xc   \\
!    9:-a2    &  C$_{s}^{4}$   &  Cc:-a2 = Bn11               &   B -2xbc  \\
!    9:-a3    &  C$_{s}^{4}$   &  Cc:-a3 = Ib11               &   I -2xb   \\
!   10:b      &  C$_{2h}^{1}$  &  P2/m:b = P12/m1             &  -P 2y     \\
!   10:c      &  C$_{2h}^{1}$  &  P2/m:c = P112/m             &  -P 2      \\
!   10:a      &  C$_{2h}^{1}$  &  P2/m:a = P2/m11             &  -P 2x     \\
!   11:b      &  C$_{2h}^{2}$  &  P21/m:b = P121/m1           &  -P 2yb    \\
!   11:c      &  C$_{2h}^{2}$  &  P21/m:c = P1121/m           &  -P 2c     \\
!   11:a      &  C$_{2h}^{2}$  &  P21/m:a = P21/m11           &  -P 2xa    \\
!   12:b1     &  C$_{2h}^{3}$  &  C2/m:b1 = C12/m1            &  -C 2y     \\
!   12:b2     &  C$_{2h}^{3}$  &  C2/m:b2 = A12/m1            &  -A 2y     \\
!   12:b3     &  C$_{2h}^{3}$  &  C2/m:b3 = I12/m1            &  -I 2y     \\
!   12:c1     &  C$_{2h}^{3}$  &  C2/m:c1 = A112/m            &  -A 2      \\
!   12:c2     &  C$_{2h}^{3}$  &  C2/m:c2 = B112/m = B2/m     &  -B 2      \\
!   12:c3     &  C$_{2h}^{3}$  &  C2/m:c3 = I112/m            &  -I 2      \\
!   12:a1     &  C$_{2h}^{3}$  &  C2/m:a1 = B2/m11            &  -B 2x     \\
!   12:a2     &  C$_{2h}^{3}$  &  C2/m:a2 = C2/m11            &  -C 2x     \\
!   12:a3     &  C$_{2h}^{3}$  &  C2/m:a3 = I2/m11            &  -I 2x     \\
!   13:b1     &  C$_{2h}^{4}$  &  P2/c:b1 = P12/c1            &  -P 2yc    \\
!   13:b2     &  C$_{2h}^{4}$  &  P2/c:b2 = P12/n1            &  -P 2yac   \\
!   13:b3     &  C$_{2h}^{4}$  &  P2/c:b3 = P12/a1            &  -P 2ya    \\
!   13:c1     &  C$_{2h}^{4}$  &  P2/c:c1 = P112/a            &  -P 2a     \\
!   13:c2     &  C$_{2h}^{4}$  &  P2/c:c2 = P112/n            &  -P 2ab    \\
!   13:c3     &  C$_{2h}^{4}$  &  P2/c:c3 = P112/b = P2/b     &  -P 2b     \\
!   13:a1     &  C$_{2h}^{4}$  &  P2/c:a1 = P2/b11            &  -P 2xb    \\
!   13:a2     &  C$_{2h}^{4}$  &  P2/c:a2 = P2/n11            &  -P 2xbc   \\
!   13:a3     &  C$_{2h}^{4}$  &  P2/c:a3 = P2/c11            &  -P 2xc    \\
!   14:b1     &  C$_{2h}^{5}$  &  P21/c:b1 = P121/c1          &  -P 2ybc   \\
!   14:b2     &  C$_{2h}^{5}$  &  P21/c:b2 = P121/n1          &  -P 2yn    \\
!   14:b3     &  C$_{2h}^{5}$  &  P21/c:b3 = P121/a1          &  -P 2yab   \\
!   14:c1     &  C$_{2h}^{5}$  &  P21/c:c1 = P1121/a          &  -P 2ac    \\
!   14:c2     &  C$_{2h}^{5}$  &  P21/c:c2 = P1121/n          &  -P 2n     \\
!   14:c3     &  C$_{2h}^{5}$  &  P21/c:c3 = P1121/b = P21/b  &  -P 2bc    \\
!   14:a1     &  C$_{2h}^{5}$  &  P21/c:a1 = P21/b11          &  -P 2xab   \\
!   14:a2     &  C$_{2h}^{5}$  &  P21/c:a2 = P21/n11          &  -P 2xn    \\
!   14:a3     &  C$_{2h}^{5}$  &  P21/c:a3 = P21/c11          &  -P 2xac   \\
!   15:b1     &  C$_{2h}^{6}$  &  C2/c:b1 = C12/c1            &  -C 2yc    \\
!   15:b2     &  C$_{2h}^{6}$  &  C2/c:b2 = A12/n1            &  -A 2yac   \\
!   15:b3     &  C$_{2h}^{6}$  &  C2/c:b3 = I12/a1            &  -I 2ya    \\
!   15:-b1    &  C$_{2h}^{6}$  &  C2/c:-b1 = A12/a1           &  -A 2ya    \\
!   15:-b2    &  C$_{2h}^{6}$  &  C2/c:-b2 = C12/n1           &  -C 2ybc   \\
!   15:-b3    &  C$_{2h}^{6}$  &  C2/c:-b3 = I12/c1           &  -I 2yc    \\
!   15:c1     &  C$_{2h}^{6}$  &  C2/c:c1 = A112/a            &  -A 2a     \\
!   15:c2     &  C$_{2h}^{6}$  &  C2/c:c2 = B112/n            &  -B 2bc    \\
!   15:c3     &  C$_{2h}^{6}$  &  C2/c:c3 = I112/b            &  -I 2b     \\
!   15:-c1    &  C$_{2h}^{6}$  &  C2/c:-c1 = B112/b = B2/b    &  -B 2b     \\
!   15:-c2    &  C$_{2h}^{6}$  &  C2/c:-c2 = A112/n           &  -A 2ac    \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!   15:-c3    &  C$_{2h}^{6}$  &  C2/c:-c3 = I112/a           &  -I 2a     \\
!   15:a1     &  C$_{2h}^{6}$  &  C2/c:a1 = B2/b11            &  -B 2xb    \\
!   15:a2     &  C$_{2h}^{6}$  &  C2/c:a2 = C2/n11            &  -C 2xbc   \\
!   15:a3     &  C$_{2h}^{6}$  &  C2/c:a3 = I2/c11            &  -I 2xc    \\
!   15:-a1    &  C$_{2h}^{6}$  &  C2/c:-a1 = C2/c11           &  -C 2xc    \\
!   15:-a2    &  C$_{2h}^{6}$  &  C2/c:-a2 = B2/n11           &  -B 2xbc   \\
!   15:-a3    &  C$_{2h}^{6}$  &  C2/c:-a3 = I2/b11           &  -I 2xb    \\
!   16        &  D$_{2}^{1}$   &  P222                &   P 2 2            \\
!   17        &  D$_{2}^{2}$   &  P2221               &   P 2c 2           \\
!   17:cab    &  D$_{2}^{2}$   &  P2122               &   P 2a 2a          \\
!   17:bca    &  D$_{2}^{2}$   &  P2212               &   P 2 2b           \\
!   18        &  D$_{2}^{3}$   &  P21212              &   P 2 2ab          \\
!   18:cab    &  D$_{2}^{3}$   &  P22121              &   P 2bc 2          \\
!   18:bca    &  D$_{2}^{3}$   &  P21221              &   P 2ac 2ac        \\
!   19        &  D$_{2}^{4}$   &  P212121             &   P 2ac 2ab        \\
!   20        &  D$_{2}^{5}$   &  C2221               &   C 2c 2           \\
!   20:cab    &  D$_{2}^{5}$   &  A2122               &   A 2a 2a          \\
!   20:bca    &  D$_{2}^{5}$   &  B2212               &   B 2 2b           \\
!   21        &  D$_{2}^{6}$   &  C222                &   C 2 2            \\
!   21:cab    &  D$_{2}^{6}$   &  A222                &   A 2 2            \\
!   21:bca    &  D$_{2}^{6}$   &  B222                &   B 2 2            \\
!   22        &  D$_{2}^{7}$   &  F222                &   F 2 2            \\
!   23        &  D$_{2}^{8}$   &  I222                &   I 2 2            \\
!   24        &  D$_{2}^{9}$   &  I212121             &   I 2b 2c          \\
!   25        &  C$_{2v}^{1}$  &  Pmm2                &   P 2 -2           \\
!   25:cab    &  C$_{2v}^{1}$  &  P2mm                &   P -2 2           \\
!   25:bca    &  C$_{2v}^{1}$  &  Pm2m                &   P -2 -2          \\
!   26        &  C$_{2v}^{2}$  &  Pmc21               &   P 2c -2          \\
!   26:ba-c   &  C$_{2v}^{2}$  &  Pcm21               &   P 2c -2c         \\
!   26:cab    &  C$_{2v}^{2}$  &  P21ma               &   P -2a 2a         \\
!   26:-cba   &  C$_{2v}^{2}$  &  P21am               &   P -2 2a          \\
!   26:bca    &  C$_{2v}^{2}$  &  Pb21m               &   P -2 -2b         \\
!   26:a-cb   &  C$_{2v}^{2}$  &  Pm21b               &   P -2b -2         \\
!   27        &  C$_{2v}^{3}$  &  Pcc2                &   P 2 -2c          \\
!   27:cab    &  C$_{2v}^{3}$  &  P2aa                &   P -2a 2          \\
!   27:bca    &  C$_{2v}^{3}$  &  Pb2b                &   P -2b -2b        \\
!   28        &  C$_{2v}^{4}$  &  Pma2                &   P 2 -2a          \\
!   28:ba-c   &  C$_{2v}^{4}$  &  Pbm2                &   P 2 -2b          \\
!   28:cab    &  C$_{2v}^{4}$  &  P2mb                &   P -2b 2          \\
!   28:-cba   &  C$_{2v}^{4}$  &  P2cm                &   P -2c 2          \\
!   28:bca    &  C$_{2v}^{4}$  &  Pc2m                &   P -2c -2c        \\
!   28:a-cb   &  C$_{2v}^{4}$  &  Pm2a                &   P -2a -2a        \\
!   29        &  C$_{2v}^{5}$  &  Pca21               &   P 2c -2ac        \\
!   29:ba-c   &  C$_{2v}^{5}$  &  Pbc21               &   P 2c -2b         \\
!   29:cab    &  C$_{2v}^{5}$  &  P21ab               &   P -2b 2a         \\
!   29:-cba   &  C$_{2v}^{5}$  &  P21ca               &   P -2ac 2a        \\
!   29:bca    &  C$_{2v}^{5}$  &  Pc21b               &   P -2bc -2c       \\
!   29:a-cb   &  C$_{2v}^{5}$  &  Pb21a               &   P -2a -2ab       \\
!   30        &  C$_{2v}^{6}$  &  Pnc2                &   P 2 -2bc         \\
!   30:ba-c   &  C$_{2v}^{6}$  &  Pcn2                &   P 2 -2ac         \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!   30:cab    &  C$_{2v}^{6}$  &  P2na                &   P -2ac 2         \\
!   30:-cba   &  C$_{2v}^{6}$  &  P2an                &   P -2ab 2         \\
!   30:bca    &  C$_{2v}^{6}$  &  Pb2n                &   P -2ab -2ab      \\
!   30:a-cb   &  C$_{2v}^{6}$  &  Pn2b                &   P -2bc -2bc      \\
!   31        &  C$_{2v}^{7}$  &  Pmn21               &   P 2ac -2         \\
!   31:ba-c   &  C$_{2v}^{7}$  &  Pnm21               &   P 2bc -2bc       \\
!   31:cab    &  C$_{2v}^{7}$  &  P21mn               &   P -2ab 2ab       \\
!   31:-cba   &  C$_{2v}^{7}$  &  P21nm               &   P -2 2ac         \\
!   31:bca    &  C$_{2v}^{7}$  &  Pn21m               &   P -2 -2bc        \\
!   31:a-cb   &  C$_{2v}^{7}$  &  Pm21n               &   P -2ab -2        \\
!   32        &  C$_{2v}^{8}$  &  Pba2                &   P 2 -2ab         \\
!   32:cab    &  C$_{2v}^{8}$  &  P2cb                &   P -2bc 2         \\
!   32:bca    &  C$_{2v}^{8}$  &  Pc2a                &   P -2ac -2ac      \\
!   33        &  C$_{2v}^{9}$  &  Pna21               &   P 2c -2n         \\
!   33:ba-c   &  C$_{2v}^{9}$  &  Pbn21               &   P 2c -2ab        \\
!   33:cab    &  C$_{2v}^{9}$  &  P21nb               &   P -2bc 2a        \\
!   33:-cba   &  C$_{2v}^{9}$  &  P21cn               &   P -2n 2a         \\
!   33:bca    &  C$_{2v}^{9}$  &  Pc21n               &   P -2n -2ac       \\
!   33:a-cb   &  C$_{2v}^{9}$  &  Pn21a               &   P -2ac -2n       \\
!   34        &  C$_{2v}^{10}$ &  Pnn2                &   P 2 -2n          \\
!   34:cab    &  C$_{2v}^{10}$ &  P2nn                &   P -2n 2          \\
!   34:bca    &  C$_{2v}^{10}$ &  Pn2n                &   P -2n -2n        \\
!   35        &  C$_{2v}^{11}$ &  Cmm2                &   C 2 -2           \\
!   35:cab    &  C$_{2v}^{11}$ &  A2mm                &   A -2 2           \\
!   35:bca    &  C$_{2v}^{11}$ &  Bm2m                &   B -2 -2          \\
!   36        &  C$_{2v}^{12}$ &  Cmc21               &   C 2c -2          \\
!   36:ba-c   &  C$_{2v}^{12}$ &  Ccm21               &   C 2c -2c         \\
!   36:cab    &  C$_{2v}^{12}$ &  A21ma               &   A -2a 2a         \\
!   36:-cba   &  C$_{2v}^{12}$ &  A21am               &   A -2 2a          \\
!   36:bca    &  C$_{2v}^{12}$ &  Bb21m               &   B -2 -2b         \\
!   36:a-cb   &  C$_{2v}^{12}$ &  Bm21b               &   B -2b -2         \\
!   37        &  C$_{2v}^{13}$ &  Ccc2                &   C 2 -2c          \\
!   37:cab    &  C$_{2v}^{13}$ &  A2aa                &   A -2a 2          \\
!   37:bca    &  C$_{2v}^{13}$ &  Bb2b                &   B -2b -2b        \\
!   38        &  C$_{2v}^{14}$ &  Amm2                &   A 2 -2           \\
!   38:ba-c   &  C$_{2v}^{14}$ &  Bmm2                &   B 2 -2           \\
!   38:cab    &  C$_{2v}^{14}$ &  B2mm                &   B -2 2           \\
!   38:-cba   &  C$_{2v}^{14}$ &  C2mm                &   C -2 2           \\
!   38:bca    &  C$_{2v}^{14}$ &  Cm2m                &   C -2 -2          \\
!   38:a-cb   &  C$_{2v}^{14}$ &  Am2m                &   A -2 -2          \\
!   39        &  C$_{2v}^{15}$ &  Abm2                &   A 2 -2c          \\
!   39:ba-c   &  C$_{2v}^{15}$ &  Bma2                &   B 2 -2c          \\
!   39:cab    &  C$_{2v}^{15}$ &  B2cm                &   B -2c 2          \\
!   39:-cba   &  C$_{2v}^{15}$ &  C2mb                &   C -2b 2          \\
!   39:bca    &  C$_{2v}^{15}$ &  Cm2a                &   C -2b -2b        \\
!   39:a-cb   &  C$_{2v}^{15}$ &  Ac2m                &   A -2c -2c        \\
!   40        &  C$_{2v}^{16}$ &  Ama2                &   A 2 -2a          \\
!   40:ba-c   &  C$_{2v}^{16}$ &  Bbm2                &   B 2 -2b          \\
!   40:cab    &  C$_{2v}^{16}$ &  B2mb                &   B -2b 2          \\
!   40:-cba   &  C$_{2v}^{16}$ &  C2cm                &   C -2c 2          \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!   40:bca    &  C$_{2v}^{16}$ &  Cc2m                &   C -2c -2c        \\
!   40:a-cb   &  C$_{2v}^{16}$ &  Am2a                &   A -2a -2a        \\
!   41        &  C$_{2v}^{17}$ &  Aba2                &   A 2 -2ac         \\
!   41:ba-c   &  C$_{2v}^{17}$ &  Bba2                &   B 2 -2bc         \\
!   41:cab    &  C$_{2v}^{17}$ &  B2cb                &   B -2bc 2         \\
!   41:-cba   &  C$_{2v}^{17}$ &  C2cb                &   C -2bc 2         \\
!   41:bca    &  C$_{2v}^{17}$ &  Cc2a                &   C -2bc -2bc      \\
!   41:a-cb   &  C$_{2v}^{17}$ &  Ac2a                &   A -2ac -2ac      \\
!   42        &  C$_{2v}^{18}$ &  Fmm2                &   F 2 -2           \\
!   42:cab    &  C$_{2v}^{18}$ &  F2mm                &   F -2 2           \\
!   42:bca    &  C$_{2v}^{18}$ &  Fm2m                &   F -2 -2          \\
!   43        &  C$_{2v}^{19}$ &  Fdd2                &   F 2 -2d          \\
!   43:cab    &  C$_{2v}^{19}$ &  F2dd                &   F -2d 2          \\
!   43:bca    &  C$_{2v}^{19}$ &  Fd2d                &   F -2d -2d        \\
!   44        &  C$_{2v}^{20}$ &  Imm2                &   I 2 -2           \\
!   44:cab    &  C$_{2v}^{20}$ &  I2mm                &   I -2 2           \\
!   44:bca    &  C$_{2v}^{20}$ &  Im2m                &   I -2 -2          \\
!   45        &  C$_{2v}^{21}$ &  Iba2                &   I 2 -2c          \\
!   45:cab    &  C$_{2v}^{21}$ &  I2cb                &   I -2a 2          \\
!   45:bca    &  C$_{2v}^{21}$ &  Ic2a                &   I -2b -2b        \\
!   46        &  C$_{2v}^{22}$ &  Ima2                &   I 2 -2a          \\
!   46:ba-c   &  C$_{2v}^{22}$ &  Ibm2                &   I 2 -2b          \\
!   46:cab    &  C$_{2v}^{22}$ &  I2mb                &   I -2b 2          \\
!   46:-cba   &  C$_{2v}^{22}$ &  I2cm                &   I -2c 2          \\
!   46:bca    &  C$_{2v}^{22}$ &  Ic2m                &   I -2c -2c        \\
!   46:a-cb   &  C$_{2v}^{22}$ &  Im2a                &   I -2a -2a        \\
!   47        &  D$_{2h}^{1}$  &  Pmmm                &  -P 2 2            \\
!   48:1      &  D$_{2h}^{2}$  &  Pnnn:1              &   P 2 2 -1n        \\
!   48:2      &  D$_{2h}^{2}$  &  Pnnn:2              &  -P 2ab 2bc        \\
!   49        &  D$_{2h}^{3}$  &  Pccm                &  -P 2 2c           \\
!   49:cab    &  D$_{2h}^{3}$  &  Pmaa                &  -P 2a 2           \\
!   49:bca    &  D$_{2h}^{3}$  &  Pbmb                &  -P 2b 2b          \\
!   50:1      &  D$_{2h}^{4}$  &  Pban:1              &   P 2 2 -1ab       \\
!   50:2      &  D$_{2h}^{4}$  &  Pban:2              &  -P 2ab 2b         \\
!   50:1cab   &  D$_{2h}^{4}$  &  Pncb:1              &   P 2 2 -1bc       \\
!   50:2cab   &  D$_{2h}^{4}$  &  Pncb:2              &  -P 2b 2bc         \\
!   50:1bca   &  D$_{2h}^{4}$  &  Pcna:1              &   P 2 2 -1ac       \\
!   50:2bca   &  D$_{2h}^{4}$  &  Pcna:2              &  -P 2a 2c          \\
!   51        &  D$_{2h}^{5}$  &  Pmma                &  -P 2a 2a          \\
!   51:ba-c   &  D$_{2h}^{5}$  &  Pmmb                &  -P 2b 2           \\
!   51:cab    &  D$_{2h}^{5}$  &  Pbmm                &  -P 2 2b           \\
!   51:-cba   &  D$_{2h}^{5}$  &  Pcmm                &  -P 2c 2c          \\
!   51:bca    &  D$_{2h}^{5}$  &  Pmcm                &  -P 2c 2           \\
!   51:a-cb   &  D$_{2h}^{5}$  &  Pmam                &  -P 2 2a           \\
!   52        &  D$_{2h}^{6}$  &  Pnna                &  -P 2a 2bc         \\
!   52:ba-c   &  D$_{2h}^{6}$  &  Pnnb                &  -P 2b 2n          \\
!   52:cab    &  D$_{2h}^{6}$  &  Pbnn                &  -P 2n 2b          \\
!   52:-cba   &  D$_{2h}^{6}$  &  Pcnn                &  -P 2ab 2c         \\
!   52:bca    &  D$_{2h}^{6}$  &  Pncn                &  -P 2ab 2n         \\
!   52:a-cb   &  D$_{2h}^{6}$  &  Pnan                &  -P 2n 2bc         \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!   53        &  D$_{2h}^{7}$  &  Pmna                &  -P 2ac 2          \\
!   53:ba-c   &  D$_{2h}^{7}$  &  Pnmb                &  -P 2bc 2bc        \\
!   53:cab    &  D$_{2h}^{7}$  &  Pbmn                &  -P 2ab 2ab        \\
!   53:-cba   &  D$_{2h}^{7}$  &  Pcnm                &  -P 2 2ac          \\
!   53:bca    &  D$_{2h}^{7}$  &  Pncm                &  -P 2 2bc          \\
!   53:a-cb   &  D$_{2h}^{7}$  &  Pman                &  -P 2ab 2          \\
!   54        &  D$_{2h}^{8}$  &  Pcca                &  -P 2a 2ac         \\
!   54:ba-c   &  D$_{2h}^{8}$  &  Pccb                &  -P 2b 2c          \\
!   54:cab    &  D$_{2h}^{8}$  &  Pbaa                &  -P 2a 2b          \\
!   54:-cba   &  D$_{2h}^{8}$  &  Pcaa                &  -P 2ac 2c         \\
!   54:bca    &  D$_{2h}^{8}$  &  Pbcb                &  -P 2bc 2b         \\
!   54:a-cb   &  D$_{2h}^{8}$  &  Pbab                &  -P 2b 2ab         \\
!   55        &  D$_{2h}^{9}$  &  Pbam                &  -P 2 2ab          \\
!   55:cab    &  D$_{2h}^{9}$  &  Pmcb                &  -P 2bc 2          \\
!   55:bca    &  D$_{2h}^{9}$  &  Pcma                &  -P 2ac 2ac        \\
!   56        &  D$_{2h}^{10}$ &  Pccn                &  -P 2ab 2ac        \\
!   56:cab    &  D$_{2h}^{10}$ &  Pnaa                &  -P 2ac 2bc        \\
!   56:bca    &  D$_{2h}^{10}$ &  Pbnb                &  -P 2bc 2ab        \\
!   57        &  D$_{2h}^{11}$ &  Pbcm                &  -P 2c 2b          \\
!   57:ba-c   &  D$_{2h}^{11}$ &  Pcam                &  -P 2c 2ac         \\
!   57:cab    &  D$_{2h}^{11}$ &  Pmca                &  -P 2ac 2a         \\
!   57:-cba   &  D$_{2h}^{11}$ &  Pmab                &  -P 2b 2a          \\
!   57:bca    &  D$_{2h}^{11}$ &  Pbma                &  -P 2a 2ab         \\
!   57:a-cb   &  D$_{2h}^{11}$ &  Pcmb                &  -P 2bc 2c         \\
!   58        &  D$_{2h}^{12}$ &  Pnnm                &  -P 2 2n           \\
!   58:cab    &  D$_{2h}^{12}$ &  Pmnn                &  -P 2n 2           \\
!   58:bca    &  D$_{2h}^{12}$ &  Pnmn                &  -P 2n 2n          \\
!   59:1      &  D$_{2h}^{13}$ &  Pmmn:1              &   P 2 2ab -1ab     \\
!   59:2      &  D$_{2h}^{13}$ &  Pmmn:2              &  -P 2ab 2a         \\
!   59:1cab   &  D$_{2h}^{13}$ &  Pnmm:1              &   P 2bc 2 -1bc     \\
!   59:2cab   &  D$_{2h}^{13}$ &  Pnmm:2              &  -P 2c 2bc         \\
!   59:1bca   &  D$_{2h}^{13}$ &  Pmnm:1              &   P 2ac 2ac -1ac   \\
!   59:2bca   &  D$_{2h}^{13}$ &  Pmnm:2              &  -P 2c 2a          \\
!   60        &  D$_{2h}^{14}$ &  Pbcn                &  -P 2n 2ab         \\
!   60:ba-c   &  D$_{2h}^{14}$ &  Pcan                &  -P 2n 2c          \\
!   60:cab    &  D$_{2h}^{14}$ &  Pnca                &  -P 2a 2n          \\
!   60:-cba   &  D$_{2h}^{14}$ &  Pnab                &  -P 2bc 2n         \\
!   60:bca    &  D$_{2h}^{14}$ &  Pbna                &  -P 2ac 2b         \\
!   60:a-cb   &  D$_{2h}^{14}$ &  Pcnb                &  -P 2b 2ac         \\
!   61        &  D$_{2h}^{15}$ &  Pbca                &  -P 2ac 2ab        \\
!   61:ba-c   &  D$_{2h}^{15}$ &  Pcab                &  -P 2bc 2ac        \\
!   62        &  D$_{2h}^{16}$ &  Pnma                &  -P 2ac 2n         \\
!   62:ba-c   &  D$_{2h}^{16}$ &  Pmnb                &  -P 2bc 2a         \\
!   62:cab    &  D$_{2h}^{16}$ &  Pbnm                &  -P 2c 2ab         \\
!   62:-cba   &  D$_{2h}^{16}$ &  Pcmn                &  -P 2n 2ac         \\
!   62:bca    &  D$_{2h}^{16}$ &  Pmcn                &  -P 2n 2a          \\
!   62:a-cb   &  D$_{2h}^{16}$ &  Pnam                &  -P 2c 2n          \\
!   63        &  D$_{2h}^{17}$ &  Cmcm                &  -C 2c 2           \\
!   63:ba-c   &  D$_{2h}^{17}$ &  Ccmm                &  -C 2c 2c          \\
!   63:cab    &  D$_{2h}^{17}$ &  Amma                &  -A 2a 2a          \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!   63:-cba   &  D$_{2h}^{17}$ &  Amam                &  -A 2 2a           \\
!   63:bca    &  D$_{2h}^{17}$ &  Bbmm                &  -B 2 2b           \\
!   63:a-cb   &  D$_{2h}^{17}$ &  Bmmb                &  -B 2b 2           \\
!   64        &  D$_{2h}^{18}$ &  Cmca                &  -C 2bc 2          \\
!   64:ba-c   &  D$_{2h}^{18}$ &  Ccmb                &  -C 2bc 2bc        \\
!   64:cab    &  D$_{2h}^{18}$ &  Abma                &  -A 2ac 2ac        \\
!   64:-cba   &  D$_{2h}^{18}$ &  Acam                &  -A 2 2ac          \\
!   64:bca    &  D$_{2h}^{18}$ &  Bbcm                &  -B 2 2bc          \\
!   64:a-cb   &  D$_{2h}^{18}$ &  Bmab                &  -B 2bc 2          \\
!   65        &  D$_{2h}^{19}$ &  Cmmm                &  -C 2 2            \\
!   65:cab    &  D$_{2h}^{19}$ &  Ammm                &  -A 2 2            \\
!   65:bca    &  D$_{2h}^{19}$ &  Bmmm                &  -B 2 2            \\
!   66        &  D$_{2h}^{20}$ &  Cccm                &  -C 2 2c           \\
!   66:cab    &  D$_{2h}^{20}$ &  Amaa                &  -A 2a 2           \\
!   66:bca    &  D$_{2h}^{20}$ &  Bbmb                &  -B 2b 2b          \\
!   67        &  D$_{2h}^{21}$ &  Cmma                &  -C 2b 2           \\
!   67:ba-c   &  D$_{2h}^{21}$ &  Cmmb                &  -C 2b 2b          \\
!   67:cab    &  D$_{2h}^{21}$ &  Abmm                &  -A 2c 2c          \\
!   67:-cba   &  D$_{2h}^{21}$ &  Acmm                &  -A 2 2c           \\
!   67:bca    &  D$_{2h}^{21}$ &  Bmcm                &  -B 2 2c           \\
!   67:a-cb   &  D$_{2h}^{21}$ &  Bmam                &  -B 2c 2           \\
!   68:1      &  D$_{2h}^{22}$ &  Ccca:1              &   C 2 2 -1bc       \\
!   68:2      &  D$_{2h}^{22}$ &  Ccca:2              &  -C 2b 2bc         \\
!   68:1ba-c  &  D$_{2h}^{22}$ &  Cccb:1              &   C 2 2 -1bc       \\
!   68:2ba-c  &  D$_{2h}^{22}$ &  Cccb:2              &  -C 2b 2c          \\
!   68:1cab   &  D$_{2h}^{22}$ &  Abaa:1              &   A 2 2 -1ac       \\
!   68:2cab   &  D$_{2h}^{22}$ &  Abaa:2              &  -A 2a 2c          \\
!   68:1-cba  &  D$_{2h}^{22}$ &  Acaa:1              &   A 2 2 -1ac       \\
!   68:2-cba  &  D$_{2h}^{22}$ &  Acaa:2              &  -A 2ac 2c         \\
!   68:1bca   &  D$_{2h}^{22}$ &  Bbcb:1              &   B 2 2 -1bc       \\
!   68:2bca   &  D$_{2h}^{22}$ &  Bbcb:2              &  -B 2bc 2b         \\
!   68:1a-cb  &  D$_{2h}^{22}$ &  Bbab:1              &   B 2 2 -1bc       \\
!   68:2a-cb  &  D$_{2h}^{22}$ &  Bbab:2              &  -B 2b 2bc         \\
!   69        &  D$_{2h}^{23}$ &  Fmmm                &  -F 2 2            \\
!   70:1      &  D$_{2h}^{24}$ &  Fddd:1              &   F 2 2 -1d        \\
!   70:2      &  D$_{2h}^{24}$ &  Fddd:2              &  -F 2uv 2vw        \\
!   71        &  D$_{2h}^{25}$ &  Immm                &  -I 2 2            \\
!   72        &  D$_{2h}^{26}$ &  Ibam                &  -I 2 2c           \\
!   72:cab    &  D$_{2h}^{26}$ &  Imcb                &  -I 2a 2           \\
!   72:bca    &  D$_{2h}^{26}$ &  Icma                &  -I 2b 2b          \\
!   73        &  D$_{2h}^{27}$ &  Ibca                &  -I 2b 2c          \\
!   73:ba-c   &  D$_{2h}^{27}$ &  Icab                &  -I 2a 2b          \\
!   74        &  D$_{2h}^{28}$ &  Imma                &  -I 2b 2           \\
!   74:ba-c   &  D$_{2h}^{28}$ &  Immb                &  -I 2a 2a          \\
!   74:cab    &  D$_{2h}^{28}$ &  Ibmm                &  -I 2c 2c          \\
!   74:-cba   &  D$_{2h}^{28}$ &  Icmm                &  -I 2 2b           \\
!   74:bca    &  D$_{2h}^{28}$ &  Imcm                &  -I 2 2a           \\
!   74:a-cb   &  D$_{2h}^{28}$ &  Imam                &  -I 2c 2           \\
!   75        &  C$_{4}^{1}$   &  P4                  &   P 4              \\
!   76        &  C$_{4}^{2}$   &  P41                 &   P 4w             \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!   77        &  C$_{4}^{3}$   &  P42                 &   P 4c             \\
!   78        &  C$_{4}^{4}$   &  P43                 &   P 4cw            \\
!   79        &  C$_{4}^{5}$   &  I4                  &   I 4              \\
!   80        &  C$_{4}^{6}$   &  I41                 &   I 4bw            \\
!   81        &  S$_{4}^{1}$   &  P-4                 &   P -4             \\
!   82        &  S$_{4}^{2}$   &  I-4                 &   I -4             \\
!   83        &  C$_{4h}^{1}$  &  P4/m                &  -P 4              \\
!   84        &  C$_{4h}^{2}$  &  P42/m               &  -P 4c             \\
!   85:1      &  C$_{4h}^{3}$  &  P4/n:1              &   P 4ab -1ab       \\
!   85:2      &  C$_{4h}^{3}$  &  P4/n:2              &  -P 4a             \\
!   86:1      &  C$_{4h}^{4}$  &  P42/n:1             &   P 4n -1n         \\
!   86:2      &  C$_{4h}^{4}$  &  P42/n:2             &  -P 4bc            \\
!   87        &  C$_{4h}^{5}$  &  I4/m                &  -I 4              \\
!   88:1      &  C$_{4h}^{6}$  &  I41/a:1             &   I 4bw -1bw       \\
!   88:2      &  C$_{4h}^{6}$  &  I41/a:2             &  -I 4ad            \\
!   89        &  D$_{4}^{1}$   &  P422                &   P 4 2            \\
!   90        &  D$_{4}^{2}$   &  P4212               &   P 4ab 2ab        \\
!   91        &  D$_{4}^{3}$   &  P4122               &   P 4w 2c          \\
!   92        &  D$_{4}^{4}$   &  P41212              &   P 4abw 2nw       \\
!   93        &  D$_{4}^{5}$   &  P4222               &   P 4c 2           \\
!   94        &  D$_{4}^{6}$   &  P42212              &   P 4n 2n          \\
!   95        &  D$_{4}^{7}$   &  P4322               &   P 4cw 2c         \\
!   96        &  D$_{4}^{8}$   &  P43212              &   P 4nw 2abw       \\
!   97        &  D$_{4}^{9}$   &  I422                &   I 4 2            \\
!   98        &  D$_{4}^{10}$  &  I4122               &   I 4bw 2bw        \\
!   99        &  C$_{4v}^{1}$  &  P4mm                &   P 4 -2           \\
!  100        &  C$_{4v}^{2}$  &  P4bm                &   P 4 -2ab         \\
!  101        &  C$_{4v}^{3}$  &  P42cm               &   P 4c -2c         \\
!  102        &  C$_{4v}^{4}$  &  P42nm               &   P 4n -2n         \\
!  103        &  C$_{4v}^{5}$  &  P4cc                &   P 4 -2c          \\
!  104        &  C$_{4v}^{6}$  &  P4nc                &   P 4 -2n          \\
!  105        &  C$_{4v}^{7}$  &  P42mc               &   P 4c -2          \\
!  106        &  C$_{4v}^{8}$  &  P42bc               &   P 4c -2ab        \\
!  107        &  C$_{4v}^{9}$  &  I4mm                &   I 4 -2           \\
!  108        &  C$_{4v}^{10}$ &  I4cm                &   I 4 -2c          \\
!  109        &  C$_{4v}^{11}$ &  I41md               &   I 4bw -2         \\
!  110        &  C$_{4v}^{12}$ &  I41cd               &   I 4bw -2c        \\
!  111        &  D$_{2d}^{1}$  &  P-42m               &   P -4 2           \\
!  112        &  D$_{2d}^{2}$  &  P-42c               &   P -4 2c          \\
!  113        &  D$_{2d}^{3}$  &  P-421m              &   P -4 2ab         \\
!  114        &  D$_{2d}^{4}$  &  P-421c              &   P -4 2n          \\
!  115        &  D$_{2d}^{5}$  &  P-4m2               &   P -4 -2          \\
!  116        &  D$_{2d}^{6}$  &  P-4c2               &   P -4 -2c         \\
!  117        &  D$_{2d}^{7}$  &  P-4b2               &   P -4 -2ab        \\
!  118        &  D$_{2d}^{8}$  &  P-4n2               &   P -4 -2n         \\
!  119        &  D$_{2d}^{9}$  &  I-4m2               &   I -4 -2          \\
!  120        &  D$_{2d}^{10}$ &  I-4c2               &   I -4 -2c         \\
!  121        &  D$_{2d}^{11}$ &  I-42m               &   I -4 2           \\
!  122        &  D$_{2d}^{12}$ &  I-42d               &   I -4 2bw         \\
!  123        &  D$_{4h}^{1}$  &  P4/mmm              &  -P 4 2            \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!  124        &  D$_{4h}^{2}$  &  P4/mcc              &  -P 4 2c           \\
!  125:1      &  D$_{4h}^{3}$  &  P4/nbm:1            &   P 4 2 -1ab       \\
!  125:2      &  D$_{4h}^{3}$  &  P4/nbm:2            &  -P 4a 2b          \\
!  126:1      &  D$_{4h}^{4}$  &  P4/nnc:1            &   P 4 2 -1n        \\
!  126:2      &  D$_{4h}^{4}$  &  P4/nnc:2            &  -P 4a 2bc         \\
!  127        &  D$_{4h}^{5}$  &  P4/mbm              &  -P 4 2ab          \\
!  128        &  D$_{4h}^{6}$  &  P4/mnc              &  -P 4 2n           \\
!  129:1      &  D$_{4h}^{7}$  &  P4/nmm:1            &   P 4ab 2ab -1ab   \\
!  129:2      &  D$_{4h}^{7}$  &  P4/nmm:2            &  -P 4a 2a          \\
!  130:1      &  D$_{4h}^{8}$  &  P4/ncc:1            &   P 4ab 2n -1ab    \\
!  130:2      &  D$_{4h}^{8}$  &  P4/ncc:2            &  -P 4a 2ac         \\
!  131        &  D$_{4h}^{9}$  &  P42/mmc             &  -P 4c 2           \\
!  132        &  D$_{4h}^{10}$ &  P42/mcm             &  -P 4c 2c          \\
!  133:1      &  D$_{4h}^{11}$ &  P42/nbc:1           &   P 4n 2c -1n      \\
!  133:2      &  D$_{4h}^{11}$ &  P42/nbc:2           &  -P 4ac 2b         \\
!  134:1      &  D$_{4h}^{12}$ &  P42/nnm:1           &   P 4n 2 -1n       \\
!  134:2      &  D$_{4h}^{12}$ &  P42/nnm:2           &  -P 4ac 2bc        \\
!  135        &  D$_{4h}^{13}$ &  P42/mbc             &  -P 4c 2ab         \\
!  136        &  D$_{4h}^{14}$ &  P42/mnm             &  -P 4n 2n          \\
!  137:1      &  D$_{4h}^{15}$ &  P42/nmc:1           &   P 4n 2n -1n      \\
!  137:2      &  D$_{4h}^{15}$ &  P42/nmc:2           &  -P 4ac 2a         \\
!  138:1      &  D$_{4h}^{16}$ &  P42/ncm:1           &   P 4n 2ab -1n     \\
!  138:2      &  D$_{4h}^{16}$ &  P42/ncm:2           &  -P 4ac 2ac        \\
!  139        &  D$_{4h}^{17}$ &  I4/mmm              &  -I 4 2            \\
!  140        &  D$_{4h}^{18}$ &  I4/mcm              &  -I 4 2c           \\
!  141:1      &  D$_{4h}^{19}$ &  I41/amd:1           &   I 4bw 2bw -1bw   \\
!  141:2      &  D$_{4h}^{19}$ &  I41/amd:2           &  -I 4bd 2          \\
!  142:1      &  D$_{4h}^{20}$ &  I41/acd:1           &   I 4bw 2aw -1bw   \\
!  142:2      &  D$_{4h}^{20}$ &  I41/acd:2           &  -I 4bd 2c         \\
!  143        &  C$_{3}^{1}$   &  P3                  &   P 3              \\
!  144        &  C$_{3}^{2}$   &  P31                 &   P 31             \\
!  145        &  C$_{3}^{3}$   &  P32                 &   P 32             \\
!  146:H      &  C$_{3}^{4}$   &  R3:H                &   R 3              \\
!  146:R      &  C$_{3}^{4}$   &  R3:R                &   P 3*             \\
!  147        &  C$_{3i}^{1}$  &  P-3                 &  -P 3              \\
!  148:H      &  C$_{3i}^{2}$  &  R-3:H               &  -R 3              \\
!  148:R      &  C$_{3i}^{2}$  &  R-3:R               &  -P 3*             \\
!  149        &  D$_{3}^{1}$   &  P312                &   P 3 2            \\
!  150        &  D$_{3}^{2}$   &  P321                &   P 3 2$''$        \\
!  151        &  D$_{3}^{3}$   &  P3112               &   P 31 2c (0 0 1)  \\
!  152        &  D$_{3}^{4}$   &  P3121               &   P 31 2$''$       \\
!  153        &  D$_{3}^{5}$   &  P3212               &   P 32 2c (0 0 -1) \\
!  154        &  D$_{3}^{6}$   &  P3221               &   P 32 2$''$       \\
!  155:H      &  D$_{3}^{7}$   &  R32:H               &   R 3 2$''$        \\
!  155:R      &  D$_{3}^{7}$   &  R32:R               &   P 3* 2           \\
!  156        &  C$_{3v}^{1}$  &  P3m1                &   P 3 -2$''$       \\
!  157        &  C$_{3v}^{2}$  &  P31m                &   P 3 -2           \\
!  158        &  C$_{3v}^{3}$  &  P3c1                &   P 3 -2$''$c      \\
!  159        &  C$_{3v}^{4}$  &  P31c                &   P 3 -2c          \\
!  160:H      &  C$_{3v}^{5}$  &  R3m:H               &   R 3 -2$''$       \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!  160:R      &  C$_{3v}^{5}$  &  R3m:R               &   P 3* -2          \\
!  161:H      &  C$_{3v}^{6}$  &  R3c:H               &   R 3 -2$''$c      \\
!  161:R      &  C$_{3v}^{6}$  &  R3c:R               &   P 3* -2n         \\
!  162        &  D$_{3d}^{1}$  &  P-31m               &  -P 3 2            \\
!  163        &  D$_{3d}^{2}$  &  P-31c               &  -P 3 2c           \\
!  164        &  D$_{3d}^{3}$  &  P-3m1               &  -P 3 2$''$        \\
!  165        &  D$_{3d}^{4}$  &  P-3c1               &  -P 3 2$''$c       \\
!  166:H      &  D$_{3d}^{5}$  &  R-3m:H              &  -R 3 2$''$        \\
!  166:R      &  D$_{3d}^{5}$  &  R-3m:R              &  -P 3* 2           \\
!  167:H      &  D$_{3d}^{6}$  &  R-3c:H              &  -R 3 2$''$c       \\
!  167:R      &  D$_{3d}^{6}$  &  R-3c:R              &  -P 3* 2n          \\
!  168        &  C$_{6}^{1}$   &  P6                  &   P 6              \\
!  169        &  C$_{6}^{2}$   &  P61                 &   P 61             \\
!  170        &  C$_{6}^{3}$   &  P65                 &   P 65             \\
!  171        &  C$_{6}^{4}$   &  P62                 &   P 62             \\
!  172        &  C$_{6}^{5}$   &  P64                 &   P 64             \\
!  173        &  C$_{6}^{6}$   &  P63                 &   P 6c             \\
!  174        &  C$_{3h}^{1}$  &  P-6                 &   P -6             \\
!  175        &  C$_{6h}^{1}$  &  P6/m                &  -P 6              \\
!  176        &  C$_{6h}^{2}$  &  P63/m               &  -P 6c             \\
!  177        &  D$_{6}^{1}$   &  P622                &   P 6 2            \\
!  178        &  D$_{6}^{2}$   &  P6122               &   P 61 2 (0 0 -1)  \\
!  179        &  D$_{6}^{3}$   &  P6522               &   P 65 2 (0 0 1)   \\
!  180        &  D$_{6}^{4}$   &  P6222               &   P 62 2c (0 0 1)  \\
!  181        &  D$_{6}^{5}$   &  P6422               &   P 64 2c (0 0 -1) \\
!  182        &  D$_{6}^{6}$   &  P6322               &   P 6c 2c          \\
!  183        &  C$_{6v}^{1}$  &  P6mm                &   P 6 -2           \\
!  184        &  C$_{6v}^{2}$  &  P6cc                &   P 6 -2c          \\
!  185        &  C$_{6v}^{3}$  &  P63cm               &   P 6c -2          \\
!  186        &  C$_{6v}^{4}$  &  P63mc               &   P 6c -2c         \\
!  187        &  D$_{3h}^{1}$  &  P-6m2               &   P -6 2           \\
!  188        &  D$_{3h}^{2}$  &  P-6c2               &   P -6c 2          \\
!  189        &  D$_{3h}^{3}$  &  P-62m               &   P -6 -2          \\
!  190        &  D$_{3h}^{4}$  &  P-62c               &   P -6c -2c        \\
!  191        &  D$_{6h}^{1}$  &  P6/mmm              &  -P 6 2            \\
!  192        &  D$_{6h}^{2}$  &  P6/mcc              &  -P 6 2c           \\
!  193        &  D$_{6h}^{3}$  &  P63/mcm             &  -P 6c 2           \\
!  194        &  D$_{6h}^{4}$  &  P63/mmc             &  -P 6c 2c          \\
!  195        &  T$_{}^{1}$    &  P23                 &   P 2 2 3          \\
!  196        &  T$_{}^{2}$    &  F23                 &   F 2 2 3          \\
!  197        &  T$_{}^{3}$    &  I23                 &   I 2 2 3          \\
!  198        &  T$_{}^{4}$    &  P213                &   P 2ac 2ab 3      \\
!  199        &  T$_{}^{5}$    &  I213                &   I 2b 2c 3        \\
!  200        &  T$_{h}^{1}$   &  Pm-3                &  -P 2 2 3          \\
!  201:1      &  T$_{h}^{2}$   &  Pn-3:1              &   P 2 2 3 -1n      \\
!  201:2      &  T$_{h}^{2}$   &  Pn-3:2              &  -P 2ab 2bc 3      \\
!  202        &  T$_{h}^{3}$   &  Fm-3                &  -F 2 2 3          \\
!  203:1      &  T$_{h}^{4}$   &  Fd-3:1              &   F 2 2 3 -1d      \\
!  203:2      &  T$_{h}^{4}$   &  Fd-3:2              &  -F 2uv 2vw 3      \\
!  204        &  T$_{h}^{5}$   &  Im-3                &  -I 2 2 3          \\
!  \hline
! \end{tabular}
! \newpage
! \begin{tabular}{|l|l|l|l|}
!  \hline
!  Number & Schoenflies & Hermann-Mauguin & Hall \\
!  \hline
!  205        &  T$_{h}^{6}$   &  Pa-3                &  -P 2ac 2ab 3      \\
!  206        &  T$_{h}^{7}$   &  Ia-3                &  -I 2b 2c 3        \\
!  207        &  O$^{1}$       &  P432                &   P 4 2 3          \\
!  208        &  O$^{2}$       &  P4232               &   P 4n 2 3         \\
!  209        &  O$^{3}$       &  F432                &   F 4 2 3          \\
!  210        &  O$^{4}$       &  F4132               &   F 4d 2 3         \\
!  211        &  O$^{5}$       &  I432                &   I 4 2 3          \\
!  212        &  O$^{6}$       &  P4332               &   P 4acd 2ab 3     \\
!  213        &  O$^{7}$       &  P4132               &   P 4bd 2ab 3      \\
!  214        &  O$^{8}$       &  I4132               &   I 4bd 2c 3       \\
!  215        &  T$_{d}^{1}$   &  P-43m               &   P -4 2 3         \\
!  216        &  T$_{d}^{2}$   &  F-43m               &   F -4 2 3         \\
!  217        &  T$_{d}^{3}$   &  I-43m               &   I -4 2 3         \\
!  218        &  T$_{d}^{4}$   &  P-43n               &   P -4n 2 3        \\
!  219        &  T$_{d}^{5}$   &  F-43c               &   F -4c 2 3        \\
!  220        &  T$_{d}^{6}$   &  I-43d               &   I -4bd 2c 3      \\
!  221        &  O$_{h}^{1}$   &  Pm-3m               &  -P 4 2 3          \\
!  222:1      &  O$_{h}^{2}$   &  Pn-3n:1             &   P 4 2 3 -1n      \\
!  222:2      &  O$_{h}^{2}$   &  Pn-3n:2             &  -P 4a 2bc 3       \\
!  223        &  O$_{h}^{3}$   &  Pm-3n               &  -P 4n 2 3         \\
!  224:1      &  O$_{h}^{4}$   &  Pn-3m:1             &   P 4n 2 3 -1n     \\
!  224:2      &  O$_{h}^{4}$   &  Pn-3m:2             &  -P 4bc 2bc 3      \\
!  225        &  O$_{h}^{5}$   &  Fm-3m               &  -F 4 2 3          \\
!  226        &  O$_{h}^{6}$   &  Fm-3c               &  -F 4c 2 3         \\
!  227:1      &  O$_{h}^{7}$   &  Fd-3m:1             &   F 4d 2 3 -1d     \\
!  227:2      &  O$_{h}^{7}$   &  Fd-3m:2             &  -F 4vw 2vw 3      \\
!  228:1      &  O$_{h}^{8}$   &  Fd-3c:1             &   F 4d 2 3 -1cd    \\
!  228:2      &  O$_{h}^{8}$   &  Fd-3c:2             &  -F 4cvw 2vw 3     \\
!  229        &  O$_{h}^{9}$   &  Im-3m               &  -I 4 2 3          \\
!  230        &  O$_{h}^{10}$  &  Ia-3d               &  -I 4bd 2c 3       \\
!  \hline
! \end{tabular}
! \end{center}
!
!EOI
