!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module modtetra
  ! Variable names taken from the GW implementation into the EXCITING code
  ! version 0.9.52 by R. Gomez-Abal.
      Implicit None
!
  !----------------------------!
  !     ordering variables     !
  !----------------------------!
  ! map from library k-point index to default k-point index
      Integer, Allocatable :: iktet2ik (:)
  ! reverse map
      Integer, Allocatable :: ik2iktet (:)
!
  !--------------------------------------!
  !     tetrahedron method variables     !
  !--------------------------------------!
  ! tetrahedron method is used for occupation numbers and Fermi energy
      Logical :: tetraocc
  ! tetrahedron method is used for optics
      Logical :: tetraopt
  ! tetrahedron method is used for dielectric function/matrix
      Logical :: tetradf
  ! use different k-point ordering that matches the one of the exciting code
      Logical :: tetrakordexc
  ! integer k-point offset
      Integer (4) :: ikloff (3)
  ! k-point offset divisor
      Integer (4) :: dkloff
  ! k-points common divisor
      Integer (4) :: dvk
  ! Number of tetrahedra
      Integer (4) :: ntet
  ! index of the k-points corresponding to the nodes of each tetrahedron
      Integer (4), Allocatable :: tnodes (:, :)
  ! weight of each tetrahedron.
      Integer (4), Allocatable :: wtet (:)
  ! volume of the tetrahedra relative to the BZ volume
      Real (8) :: tvol
  ! parameter specifying smalles diagonal in generic tetrahedron
      Integer (4) :: mnd
!
  !---------------------------------!
  !     q-dependent convolution     !
  !---------------------------------!
  ! number of the tetrahedra linked to by the corresponding q vector
      Integer (4), Allocatable :: link (:), kqid (:, :)
  ! q-points common divisor
      Integer (4) dvq
!
Contains
!
!BOP
! !ROUTINE: rtorat
! !INTERFACE:
      Subroutine rtorat (eps, n, x, k, div)
! !DESCRIPTION:
!   This subroutine factorizes the real coordinates of a vector {\bf x}.
!   The output is an integer vector {\bf k}, such that $k(i)/{\rm div}=x(i)$
!   and
!   $$ |x(i)-k(i)/{\rm div}| < {\rm eps} $$.
!
! !REVISION HISTORY:
!   Created July 2008 by Sagmeister
!EOP
!BOC
         Implicit None
    ! arguments
         Real (8), Intent (In) :: eps
         Integer (4), Intent (In) :: n
         Real (8), Intent (In) :: x (n)
         Integer (4), Intent (Out) :: div
         Integer (4), Intent (Out) :: k (n)
    ! local variables
         Integer :: maxint
         Real (8) :: dx
         maxint = Nint (1.d0/eps) / 10
         Do div = 1, maxint
            k (:) = Nint (dble(div)*x(:))
            dx = maxval (Abs(dble(k)/dble(div)-x))
            If (dx .Lt. eps) Exit
         End Do
         If (dx .Ge. eps) Then
            Write (*,*)
            Write (*, '("Error(modtetra:rtorat): factorization failed")&
           &')
            Write (*, '(" maximum integer :",i12)') maxint
            Write (*, '(" tolerance       :",g18.10)') eps
            Write (*, '(" deviation       :",g18.10)') dx
            Write (*,*)
            Stop
         End If
         If (dx .Gt. 1.d-12) Then
            Write (*,*)
            Write (*, '("Warning(modtetra:rtorat): small deviation in f&
           &actorization")')
            Write (*, '(" maximum deviation :",g18.10)') dx
            Write (*,*)
         End If
      End Subroutine rtorat
!EOC
!
!
      Subroutine r3fraction (r, n, d)
         Implicit None
    ! arguments
         Real (8), Intent (In) :: r (3)
         Integer, Intent (Out) :: n (3), d
    ! parameters
         Real (8), Parameter :: eps = 1.d-5, eps2 = 1.d-3
         Call rtorat (eps, 3, r, n, d)
    ! check factorization
         If ((sum(Abs(r)) .Lt. eps2) .And. (sum(Abs(r)) .Gt. 0.d0)) &
        & Then
            Write (*,*)
            Write (*, '("Warning(modtetra:r3fraction): very small offse&
           &t:")')
            Write (*, '(" kgen and related routines might fail")')
            Write (*,*)
         End If
      End Subroutine r3fraction
!
!
      Subroutine geniktetmap (eps, nppt, ngridp, vploff, vpllib, vpl, &
     & ipmap)
         Implicit None
    ! arguments
         Real (8), Intent (In) :: eps
         Integer, Intent (In) :: nppt
         Integer, Intent (In) :: ngridp (3)
         Real (8), Intent (In) :: vploff (3)
         Real (8), Intent (In) :: vpl (3, nppt)
         Real (8), Intent (In) :: vpllib (3, nppt)
         Integer, Intent (In) :: ipmap (0:ngridp(1)-1, 0:ngridp(2)-1, &
        & 0:ngridp(3)-1)
    ! local variables
         Integer :: ip, ipd, iv (3)
         If (allocated(iktet2ik)) deallocate (iktet2ik)
         Allocate (iktet2ik(nppt))
         If (allocated(ik2iktet)) deallocate (ik2iktet)
         Allocate (ik2iktet(nppt))
         Do ip = 1, nppt
       ! grid coordinates of library k-point
            iv (:) = Nint (vpllib(:, ip)*ngridp-vploff(:))
       ! index in default p-point set
            ipd = ipmap (iv(1), iv(2), iv(3))
       ! map from library to default
            iktet2ik (ip) = ipd
       ! reverse map
            ik2iktet (ipd) = ip
       ! check maps
            If (sum(Abs(vpl(:, ipd)-vpllib(:, ip))) .Ge. eps) Then
               Write (*,*)
               Write (*, '("Error(modtetra:geniktetmap): k-point mappin&
              &g between")')
               Write (*, '(" set of library and default set failed:")')
               Write (*, '(" library k-point       :",3g18.10)') vpllib &
              & (:, ip)
               Write (*, '(" mapped default k-point:",3g18.10)') vpl &
              & (:, ipd)
               Write (*, '(" map from library to default:",2i8)') ip, &
              & ipd
               Write (*,*)
               Stop
            End If
         End Do
      End Subroutine geniktetmap
!
!
      Subroutine writeiktetmap (filext, nkpt)
         Implicit None
    ! arguments
         Character (*), Intent (In) :: filext
         Integer, Intent (In) :: nkpt
    ! local variables
         Integer :: ik, un
         un = 771
         Open (un, File='KTETMAP'//trim(filext), Form='formatted', &
        & Action='write', Status='replace')
         Write (un, '(i9,a)') nkpt, ' : nkpt; k-point, iktet2ik below'
         Do ik = 1, nkpt
            Write (un, '(2i9)') ik, iktet2ik (ik)
         End Do
         Close (un)
         Open (un, File='KTETIMAP'//trim(filext), Form='formatted', &
        & Action='write', Status='replace')
         Write (un, '(i9,a)') nkpt, ' : nkpt; k-point, ik2iktet below'
         Do ik = 1, nkpt
            Write (un, '(2i9)') ik, ik2iktet (ik)
         End Do
         Close (un)
      End Subroutine writeiktetmap
!
!
      Subroutine genkpts_tet (filext, eps, bvec, maxsymcrys, nsymcrys, &
     & lsplsymc, symlat, reducek, ngridk, vkloff, nkpt, ikmap, vkl, &
     & wkpt)
         Implicit None
    ! arguments
         Character (*), Intent (In) :: filext
         Real (8), Intent (In) :: eps
         Real (8), Intent (In) :: bvec (3, 3)
         Integer, Intent (In) :: maxsymcrys
         Integer, Intent (In) :: nsymcrys
         Integer, Intent (In) :: lsplsymc (maxsymcrys)
         Integer, Intent (In) :: symlat (3, 3, 48)
         Logical, Intent (In) :: reducek
         Integer, Intent (In) :: ngridk (3)
         Real (8), Intent (In) :: vkloff (3)
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: ikmap (0:ngridk(1)-1, 0:ngridk(2)-1, &
        & 0:ngridk(3)-1)
         Real (8), Intent (In) :: vkl (3, &
        & ngridk(1)*ngridk(2)*ngridk(3))
         Real (8), Intent (In) :: wkpt (ngridk(1)*ngridk(2)*ngridk(3))
    ! local variables
         Real (8), Parameter :: epsvkloff = 1.d-5
         Integer :: isym, lspl, i1, i2, nsymcryst
         Integer :: ik, ikd, nkptlib
         Real (8) :: wkptlib
         Integer, Allocatable :: symc (:, :, :)
         Integer, Allocatable :: ivk (:, :)
         Integer, Allocatable :: indirkp (:)
         Integer, Allocatable :: iwkp (:)
         Real (8), Allocatable :: vkllib (:, :)
         If (nsymcrys .Gt. 48) Then
            Write (*,*)
            Write (*, '("Error(modtetra:genkpts_tet): number of crystal&
           & symmetries > 48:")')
            Write (*, '(" This does not work with the k-point generatio&
           &n of")')
            Write (*, '(" the linear tetrahedron method.")')
            Write (*,*)
            Stop
         End If
    ! switch to exciting interface
         If (tetrakordexc) Call tetrasetifc ('exciting')
    ! suppress debug output in tetrahedron integration library (0)
         Call tetrasetdbglv (0)
    ! safer pointer handling in tetrahedron integration library (1)
         Call tetrasetpointerhandling (1)
    ! set resonance type (1...resonant weights)
         Call tetrasetresptype (1)
    ! set treatment of q-shifted k-mesh
         Call tetrasetkplusq (.True.)
    ! report interface parameters
         Call tetrareportsettings
    ! generate fraction for k-point offset
         Call rtorat (epsvkloff, 3, vkloff, ikloff, dkloff)
    ! get rotational part of crystal symmetries
         Allocate (symc(3, 3, nsymcrys))
         Do isym = 1, nsymcrys
            lspl = lsplsymc (isym)
       ! transpose of rotation for use with the library
            Do i1 = 1, 3
               Do i2 = 1, 3
                  symc (i1, i2, isym) = symlat (i2, i1, lspl)
               End Do
            End Do
         End Do
         nsymcryst = 1
         If (reducek) nsymcryst = nsymcrys
    ! allocate local variables
         Allocate (indirkp(ngridk(1)*ngridk(2)*ngridk(3)))
         Allocate (iwkp(ngridk(1)*ngridk(2)*ngridk(3)))
         Allocate (ivk(3, ngridk(1)*ngridk(2)*ngridk(3)))
         Allocate (vkllib(3, ngridk(1)*ngridk(2)*ngridk(3)))
    ! allocate weights
         If (allocated(wtet)) deallocate (wtet)
         Allocate (wtet(1:ngridk(1)*ngridk(2)*ngridk(3)*6))
         wtet (:) = 0
    ! allocate nodes
         If (allocated(tnodes)) deallocate (tnodes)
         Allocate (tnodes(1:4, 1:ngridk(1)*ngridk(2)*ngridk(3)*6))
         tnodes (:, :) = 0
    ! generate k-point set by using library-routine
         Call kgen (bvec, nsymcryst, symc, ngridk, ikloff, dkloff, &
        & nkptlib, ivk, dvk, indirkp, iwkp, ntet, tnodes, wtet, tvol, &
        & mnd)
         If (nkptlib .Ne. nkpt) Then
            Write (*,*)
            Write (*, '("Error(modtetra:genkpts_tet): k-point set incon&
           &sistency for tetrahedron method")')
            Write (*, '(" differring number of k-points (library/default)",2i8)') nkptlib, nkpt
            Write (*,*)
            Stop
         End If
    ! k-point in lattice coordinates
         Do ik = 1, nkpt
            vkllib (:, ik) = dble (ivk(:, ik)) / dble (dvk)
         End Do
    ! generate map between the k-point set of the library and the default one
         Call geniktetmap (eps, nkpt, ngridk, vkloff, vkllib, vkl, &
        & ikmap)
    ! check weights of k-points
         Do ik = 1, nkpt
            ikd = iktet2ik (ik)
            wkptlib = iwkp (ik) / dble (ngridk(1)*ngridk(2)*ngridk(3))
            If (Abs(wkpt(ikd)-wkptlib) .Gt. eps) Then
               Write (*,*)
               Write (*, '("Error(modtetra:genkpts_tet): differring wei&
              &ghts:")')
               Write (*, '(" k-point (default)",i8,3g18.10)') ikd, vkl &
              & (:, ikd)
               Write (*, '(" weight (default)",g18.10)') wkpt (ikd)
               Write (*, '(" weight (library)",g18.10)') wkptlib
               Write (*,*)
               Stop
            End If
         End Do
         Deallocate (indirkp, iwkp, ivk, vkllib)
    ! write maps
         Call writeiktetmap (filext, nkpt)
      End Subroutine genkpts_tet
!
!
!BOP
! !ROUTINE: gentetlink
! !INTERFACE:
      Subroutine gentetlink (vpl, tqw, eps, bvec, ngridk, vkloff, nkpt, &
     & nkptnr, vklnr, ikmapnr)
! !DESCRIPTION:
!   Generates an array connecting the tetrahedra of the $\mathbf{k}$-point with
!   the ones of the  $\mathbf{k}+\mathbf{q}$-point. Interface routine
!   referencing the {\tt libbzint} library of Ricardo Gomez-Abal.
!
! !REVISION HISTORY:
!   Created January 2008 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Real (8), Intent (In) :: vpl (3)
         Integer, Intent (In) :: tqw
         Real (8), Intent (In) :: eps
         Real (8), Intent (In) :: bvec (3, 3)
         Integer, Intent (In) :: ngridk (3)
         Real (8), Intent (In) :: vkloff (3)
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: nkptnr
         Real (8), Intent (In) :: vklnr (3, &
        & ngridk(1)*ngridk(2)*ngridk(3))
         Integer, Intent (In) :: ikmapnr (0:ngridk(1)-1, 0:ngridk(2)-1, &
        & 0:ngridk(3)-1)
    ! local variables
         Real (8), Parameter :: epscomm = 1.d-5
         Real (8) :: vr (3)
         Integer :: j, iv (3), iqnr
         Logical :: tqg
         Integer, Allocatable :: ivkt (:, :), ivqt (:, :), tnodest (:, &
        & :), wtett (:)
         Real (8), External :: r3taxi
         tqg = sum (Abs(vpl)) .Lt. eps
    ! get index to reducible q-point which is commensurate to k-point set
         vr (:) = vpl (:) * ngridk (:)
         Call r3frac (eps, vr, iv)
         If (sum(Abs(vr/ngridk)) .Gt. epscomm) Then
            Write (*,*)
            Write (*, '("Error(gentetlink): q-point not commensurate wi&
           &th k-point set")')
            Write (*, '(" which is required for tetrahedron method")')
            Write (*, '(" commensurability tolerance: ",g18.10)') &
           & epscomm
            Write (*, '(" q-point (latt. coords.)   : ",3g18.10)') vpl
            Write (*, '(" deviation                 : ",3g18.10)') vr / &
           & ngridk (:)
            Write (*, '(" minimum nonzero coords.   : ",3g18.10)') 1.d0 &
           & / ngridk (:)
            Write (*,*)
            Call terminate
         End If
         iqnr = ikmapnr (iv(1), iv(2), iv(3))
    ! cross check q-point again
         vr (:) = vklnr (:, iqnr) - vkloff (:) / ngridk (:)
         If (Abs(r3taxi(vpl, vr)) .Gt. epscomm) Then
            Write (*,*)
            Write (*, '("Error(gentetlink): specified q-point does not &
           &match derived q-point on grid")')
            Write (*, '(" specified q-point :",3g18.10)') vpl
            Write (*, '(" derived q-point   :",3g18.10)') vr
            Write (*, '(" non-reduced index :",i6)') iqnr
            Write (*,*)
            Call terminate
         End If
         Write (*,*)
         Write (*, '("Info(gentetlink): q-point on grid")')
         Write (*, '(" q-point           :",3g18.10)') vr
         Write (*, '(" non-reduced index :",i6)') iqnr
         Write (*,*)
    ! check if k-point set is not reduced for q-point different from Gamma point
         If ((nkpt .Ne. nkptnr) .And. ( .Not. tqg)) Then
            Write (*,*)
            Write (*, '("Error(gentetlink): k-point set is reduced by s&
           &ymmetries and q-point is not Gamma point")')
            Write (*,*)
            Call terminate
         End If
    ! allocate link array
         If (allocated(link)) deallocate (link)
         Allocate (link(6*nkptnr))
    ! quick return for Gamma q-point
         If (tqg) Then
            Forall (j=1:6*nkptnr) link (j) = j
            Return
         End If
    ! allocate local arrays
         Allocate (ivkt(3, nkptnr), ivqt(3, nkptnr))
         Allocate (wtett(6*nkptnr), tnodest(4, 6*nkptnr))
    ! generate fraction for k-point offset
         Call r3fraction (vkloff, ikloff, dkloff)
    ! call to library routine (generate link array, nodes and weights)
         Call kqgen_exciting (bvec, ngridk, ikloff, dkloff, nkpt, iqnr, &
        & ivkt, ivqt, dvk, dvq, ntet, tnodest, wtett, link, tvol)
         If (tqw .Ne. 0) Then
       ! use weights and nodes from kqgen-routine
            tnodes (:, :) = tnodest (:, :)
            wtet (:) = wtett (:)
         Else If (sum(Abs(vpl)) .Gt. eps) Then
            Write (*,*)
            Write (*, '("Warning(gentetlink): using WTET and TNODES arr&
           &ays from KGEN routine")')
            Write (*, '(" but arrays from KQGEN_EXCITING are differring&
           & for")')
            Write (*, '(" non-Gamma q-point: ",3g18.10)') vpl
            Write (*,*)
         End If
    ! deallocate local arrays
         Deallocate (ivkt, ivqt)
         Deallocate (wtett, tnodest)
      End Subroutine gentetlink
!EOC
!
!
      Subroutine fermitetifc (nkpt, nst, eval, chgval, spinpol, efermi, &
     & fermidos)
         Implicit None
    ! arguments
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: nst
         Real (8), Intent (In) :: eval (nst, nkpt)
         Real (8), Intent (In) :: chgval
         Logical, Intent (In) :: spinpol
         Real (8), Intent (Out) :: efermi
         Real (8), Intent (Out) :: fermidos
    ! local variables
         Integer :: ik
         Real (8), Allocatable :: evallib (:, :)
         Allocate (evallib(nst, nkpt))
    ! reorder energies to library order
         Do ik = 1, nkpt
            evallib (:, ik) = eval (:, iktet2ik(ik))
         End Do
    ! call to library routine
         Call fermitet (nkpt, nst, evallib, ntet, tnodes, wtet, tvol, &
        & chgval, spinpol, efermi, fermidos, .False.)
         Deallocate (evallib)
      End Subroutine fermitetifc
!
!
      Subroutine tetiwifc (nkpt, nst, eval, efermi, occ)
         Implicit None
    ! arguments
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: nst
         Real (8), Intent (In) :: eval (nst, nkpt)
         Real (8), Intent (In) :: efermi
         Real (8), Intent (Out) :: occ (nst, nkpt)
    ! local variables
         Integer :: ik, ikd
         Real (8), Allocatable :: evallib (:, :), occt (:)
         Allocate (evallib(nst, nkpt), occt(nst))
    ! reorder energies to library order
         Do ik = 1, nkpt
            evallib (:, ik) = eval (:, iktet2ik(ik))
         End Do
    ! call to library routine
         Call tetiw (nkpt, ntet, nst, evallib, tnodes, wtet, tvol, &
        & efermi, occ)
    ! reorder occupation numbers to default order
         Do ik = 1, nkpt
            ikd = iktet2ik (ik)
            occt (:) = occ (:, ikd)
            occ (:, ikd) = occ (:, ik)
            occ (:, ik) = occt (:)
         End Do
         Deallocate (evallib, occt)
      End Subroutine tetiwifc
!
!
      Subroutine tetcwifc (nkpt, nst, eval, efermi, w, ifreq, cw)
         Implicit None
    ! arguments
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: nst
         Real (8), Intent (In) :: eval (nst, nkpt)
         Real (8), Intent (In) :: efermi
         Real (8), Intent (In) :: w
         Integer, Intent (In) :: ifreq
         Real (8), Intent (Out) :: cw (nst, nst, nkpt)
    ! local variables
         Integer :: ik, ikd
         Real (8), Allocatable :: evallib (:, :), cwt (:, :)
         Allocate (evallib(nst, nkpt), cwt(nst, nst))
    ! reorder energies to library order
         Do ik = 1, nkpt
            evallib (:, ik) = eval (:, iktet2ik(ik))
         End Do
    ! call to library routine
         Call tetcw (nkpt, ntet, nst, wtet, evallib, tnodes, link, &
        & tvol, efermi, w, ifreq, cw)
    ! reorder convolution weights to default order
         Do ik = 1, nkpt
            ikd = iktet2ik (ik)
            cwt (:, :) = cw (:, :, ikd)
            cw (:, :, ikd) = cw (:, :, ik)
            cw (:, :, ik) = cwt (:, :)
         End Do
         Deallocate (evallib, cwt)
      End Subroutine tetcwifc
!
!
      Subroutine tetcwifc_1k (ik, nkpt, nst, eval, efermi, w, ifreq, &
     & cw)
         Implicit None
    ! arguments
         Integer, Intent (In) :: ik
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: nst
         Real (8), Intent (In) :: eval (nst, nkpt)
         Real (8), Intent (In) :: efermi
         Real (8), Intent (In) :: w
         Integer, Intent (In) :: ifreq
         Real (8), Intent (Out) :: cw (nst, nst)
    ! local variables
         Integer :: ikt, iklib
         Real (8), Allocatable :: evallib (:, :)
         Allocate (evallib(nst, nkpt))
    ! reorder energies to library order
         Do ikt = 1, nkpt
            evallib (:, ikt) = eval (:, iktet2ik(ikt))
         End Do
    ! call to library routine
         iklib = ik2iktet (ik)
         Call tetcw_1k (iklib, nkpt, ntet, nst, wtet, evallib, tnodes, &
        & link, tvol, efermi, w, ifreq, cw)
         Deallocate (evallib)
      End Subroutine tetcwifc_1k
!
!
      Subroutine tetcwifc_1kbc (ik, ist, jst, nkpt, nst, eval, efermi, &
     & w, ifreq, cw)
         Implicit None
    ! arguments
         Integer, Intent (In) :: ik, ist, jst
         Integer, Intent (In) :: nkpt
         Integer, Intent (In) :: nst
         Real (8), Intent (In) :: eval (nst, nkpt)
         Real (8), Intent (In) :: efermi
         Real (8), Intent (In) :: w
         Integer, Intent (In) :: ifreq
         Real (8), Intent (Out) :: cw
    ! local variables
         Integer :: ikt, iklib
         Real (8), Allocatable :: evallib (:, :)
         Allocate (evallib(nst, nkpt))
    ! reorder energies to library order
         Do ikt = 1, nkpt
            evallib (:, ikt) = eval (:, iktet2ik(ikt))
         End Do
    ! call to library routine
         iklib = ik2iktet (ik)
         Call tetcw_1kbc (iklib, ist, jst, nkpt, ntet, nst, wtet, &
        & evallib, tnodes, link, tvol, efermi, w, ifreq, cw)
         Deallocate (evallib)
      End Subroutine tetcwifc_1kbc
!
End Module modtetra
