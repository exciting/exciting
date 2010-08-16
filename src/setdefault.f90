
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:
Subroutine setdefault
! !USES:
      Use modinput
      Use modmain
#ifdef TETRA
      Use modtetra
#endif
#ifdef XS
      Use modmpi, Only: rank
      Use modxs
#endif
      Use sclcontroll
! !DESCRIPTION:
!   Sets default values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Additional parmeters for excited states and tetrahedron method
!     2004-2008 (Sagmeister)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, js, ia, ja, ias
      Integer :: i, l, iv, iostat
      Character (256) :: bname
!------------------------!
!     default values     !
!------------------------!
      ntasks = 1
      tasks (1) = - 1
      ngridq (:) = 1
      nspecies = 0
      epsedirac = 1.d-11
      epspotatom = 1.d-6
      bflmt (:, :, :) = 0.d0
      scrpath = './'
      nvp1d = 2
      tarpack = .False.
      tlapack = .True.
      tdiis = .False.
      tjdqz = .False.
      diisfirstscl = 10
      lowesteval = - 1.d0
      packedmatrixstorage = .True.
      epsarpack = 1e-8
      maxncv = 200
      If (allocated(vvlp1d)) deallocate (vvlp1d)
      Allocate (vvlp1d(3, nvp1d))
      vvlp1d (:, 1) = 0.d0
      vvlp1d (:, 2) = 1.d0
      npp1d = 200
      vclp2d (:, :) = 0.d0
      vclp2d (1, 2) = 1.d0
      vclp2d (2, 3) = 1.d0
      np2d (:) = 40
      vclp3d (:, :) = 0.d0
      vclp3d (1, 2) = 1.d0
      vclp3d (2, 3) = 1.d0
      vclp3d (3, 4) = 1.d0
      np3d (:) = 20
      wdos (1) = - 0.5d0
      wdos (2) = 0.5d0
      nphwrt = 1
      If (allocated(vqlwrt)) deallocate (vqlwrt)
      Allocate (vqlwrt(3, nphwrt))
      vqlwrt (:, :) = 0.d0
      notelns = 0
      nkstlist = 1
      kstlist (:, 1) = 1
      ldapu = 0
      llu (:) = - 1
      ujlu (:, :) = 0.d0
      tseqit = .False.
      nseqit = 40
      tauseq = 0.1d0
#ifdef TETRA
! tetrahedron method variables
      tetraocc = .False.
      tetraopt = .False.
      tetradf = .False.
      tetrakordexc = .False.
#endif
#ifdef XS
      nbfce = - 1
      nafce = - 1
      nbfbse = - 1
      nafbse = - 1
#endif
End Subroutine
!EOC
