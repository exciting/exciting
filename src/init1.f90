!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: init1
! !INTERFACE:
!
!
Subroutine init1
! !USES:
      Use modinput
      Use modmain
      Use modmpi, only: terminate


#ifdef TETRA
      Use modtetra
#endif
#ifdef XS
      Use modxs
#endif
      use modfvsystem
      use sirius_init, only: use_sirius_gkvec
      use sirius_api,  only: setup_sirius_gs_handler, get_gkvec_arrays_sirius, &
                             get_max_num_gkvec_sirius

! !DESCRIPTION:
!   Generates the $k$-point set and then allocates and initialises global
!   variables which depend on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, is, ia, ias, io, ilo
      Integer :: i1, i2, i3, ispn, iv (3)
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Integer :: n1, n2, n3, nonzcount
      Integer :: i, maxi
      Real (8) :: vl (3), vc (3), boxl (3, 4), lambda
      Real (8) :: ts0, ts1
      real (8) :: vl1(3),vl2(3),vc1(3),vc2(3)
      real (8) :: d1,d2,d12,t1,t2
      Real (8) :: blen(3), lambdab
      logical :: wannierband
      ! +/-1 for sign of spin-dependent term
      real (8) :: sign

! external functions
      Complex (8) gauntyry
      External gauntyry
!
      call stopwatch("exciting:init1", 1)

      wannierband = .false.
      if( associated( input%properties)) then
        if( associated( input%properties%bandstructure)) wannierband = input%properties%bandstructure%wannier
      end if
      Call timesec (ts0)
!
!---------------------!
!     k-point set     !
!---------------------!
!
      If (task.Eq.100) Then
!
!     3D fermisurface plot
!
        if (associated(input%properties%fermisurfaceplot%plot3d)) then
          np3d(:) = input%properties%fermisurfaceplot%plot3d%box%grid(:)
          vclp3d(:,1) = input%properties%fermisurfaceplot%plot3d%box%origin%coord(:)
          vclp3d(:,2) = input%properties%fermisurfaceplot%plot3d%box%pointarray(1)%point%coord(:)
          vclp3d(:,3) = input%properties%fermisurfaceplot%plot3d%box%pointarray(2)%point%coord(:)
          vclp3d(:,4) = input%properties%fermisurfaceplot%plot3d%box%pointarray(3)%point%coord(:)
        else
          np3d(:)=(/20,20,20/)
          vclp3d(:,1)=(/0.d0,0.d0,0.d0/)
          vclp3d(:,2)=(/1.d0,0.d0,0.d0/)
          vclp3d(:,3)=(/0.d0,1.d0,0.d0/)
          vclp3d(:,4)=(/0.d0,0.d0,1.d0/)
        end if
        input%groundstate%ngridk(:)=np3d(:)
      End If

      If (((task .Eq. 20) .Or. (task .Eq. 21)).and..not.(wannierband)) Then
! for band structure plots generate k-points along a line
         nvp1d = size &
        & (input%properties%bandstructure%plot1d%path%pointarray)
         npp1d = input%properties%bandstructure%plot1d%path%steps
         If (allocated(dvp1d)) deallocate (dvp1d)
         Allocate (dvp1d(nvp1d))
         If (allocated(vplp1d)) deallocate (vplp1d)
         Allocate (vplp1d(3, npp1d))
         If (allocated(dpp1d)) deallocate (dpp1d)
         Allocate (dpp1d(npp1d))
         Call connect (bvec, nvp1d, npp1d, input%properties%bandstructure%plot1d%path%pointarray, &
           vplp1d, dvp1d, dpp1d)
!
         nkpt = input%properties%bandstructure%plot1d%path%steps
         nkpt_ptr => nkpt
         If (allocated(vkl)) deallocate (vkl)
         Allocate (vkl(3, nkpt))
         vkl_ptr => vkl
         If (allocated(vkc)) deallocate (vkc)
         Allocate (vkc(3, nkpt))
         Do ik = 1, nkpt
            vkl (:, ik) = vplp1d (:, ik)
            Call r3mv (bvec, vkl(:, ik), vkc(:, ik))
         End Do
      Else If (task .Eq. 25) Then
! effective mass calculation
         nkpt = (2*input%properties%masstensor%ndspem+1) ** 3
         nkpt_ptr => nkpt
         If (allocated(ivk)) deallocate (ivk)
         Allocate (ivk(3, nkpt))
         If (allocated(vkl)) deallocate (vkl)
         Allocate (vkl(3, nkpt))
         vkl_ptr => vkl
         If (allocated(vkc)) deallocate (vkc)
         Allocate (vkc(3, nkpt))
! map vector to [0,1)
         Call r3frac (input%structure%epslat, &
        & input%properties%masstensor%vklem, iv)
         ik = 0
         Do i3 = - input%properties%masstensor%ndspem, &
        & input%properties%masstensor%ndspem
            Do i2 = - input%properties%masstensor%ndspem, &
           & input%properties%masstensor%ndspem
               Do i1 = - input%properties%masstensor%ndspem, &
              & input%properties%masstensor%ndspem
                  ik = ik + 1
                  ivk (1, ik) = i1
                  ivk (2, ik) = i2
                  ivk (3, ik) = i3
                  vc (1) = dble (i1)
                  vc (2) = dble (i2)
                  vc (3) = dble (i3)
                  vc (:) = vc (:) * input%properties%masstensor%deltaem
                  Call r3mv (binv, vc, vl)
                  vkl (:, ik) = input%properties%masstensor%vklem(:) + &
                 & vl (:)
                  Call r3mv (bvec, vkl(:, ik), vkc(:, ik))
               End Do
            End Do
         End Do
!
      Else If (task .Eq. 42) Then
         !
         !          2D spin texture plot
         !
           np2d        = input%properties%spintext%plot2d%parallelogram%grid
           vclp2d(:,1) = input%properties%spintext%plot2d%parallelogram%origin%coord(:)
           vclp2d(:,2) = input%properties%spintext%plot2d%parallelogram%pointarray(1)%point%coord(:)
           vclp2d(:,3) = input%properties%spintext%plot2d%parallelogram%pointarray(2)%point%coord(:)
            ! generate 2D grid of k-points
           vl1(:) = vclp2d(:,2)-vclp2d(:,1)
           vl2(:) = vclp2d(:,3)-vclp2d(:,1)
           vc1(:) = matmul(bvec,vl1)
           vc2(:) = matmul(bvec,vl2)
           d1  = sqrt(dot_product(vc1,vc1))
           d2  = sqrt(dot_product(vc2,vc2))
           d12 = dot_product(vc1,vc2)/(d1*d2)
           if ( (d1.lt.input%structure%epslat) .or. (d2.lt.input%structure%epslat) ) then
             write (*,*)
             write (*, '("Error(fermisurf): zero length plotting vectors")')
             write (*,*)
             call terminate
           end if
           if ((1.d0-d12 <= input%structure%epslat) .or. (1.d0+d12 <= input%structure%epslat)) then
              write (*,*)
              write (*, '("Error(fermisurf): zero angle between vectors defining the parallelogram")')
              write (*,*)
              call terminate
           end if
           nkpt = product(np2d)
           If (allocated(vkl)) deallocate (vkl)
           Allocate (vkl(3, nkpt))
           If (allocated(vkc)) deallocate (vkc)
           Allocate (vkc(3, nkpt))
           ik = 0
           Do i1 = 0, np2d(1)-1
              Do i2 = 0, np2d(2)-1
                 ik = ik+1
                 t1 = dble(i1)/dble(np2d(1))
                 t2 = dble(i2)/dble(np2d(2))
                 vkl(:,ik) = t1*vl1(:)+t2*vl2(:)+vclp2d(:,1)
                 Call r3mv(bvec,vkl(:, ik),vkc(:, ik))
              End Do
           End Do
      Else If (task .Eq. 101) Then
!
!          2D fermisurface plot
!
           np2d(:) = input%properties%fermisurfaceplot%plot2d%parallelogram%grid(:)
           vclp2d(:,1) = input%properties%fermisurfaceplot%plot2d%parallelogram%origin%coord(:)
           vclp2d(:,2) = input%properties%fermisurfaceplot%plot2d%parallelogram%pointarray(1)%point%coord(:)
           vclp2d(:,3) = input%properties%fermisurfaceplot%plot2d%parallelogram%pointarray(2)%point%coord(:)
           ! generate 2D grid of k-points
           vl1(:) = vclp2d(:,2)-vclp2d(:,1)
           vl2(:) = vclp2d(:,3)-vclp2d(:,1)
           vc1(:) = bvec(:,1)*vl1(1)+bvec(:,2)*vl1(2)+bvec(:,3)*vl1(3)
           vc2(:) = bvec(:,1)*vl2(1)+bvec(:,2)*vl2(2)+bvec(:,3)*vl2(3)
           d1 = sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
           d2 = sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
           d12 = (vc1(1)*vc2(1)+vc1(2)*vc2(2)+vc1(3)*vc2(3))/(d1*d2)
           if ( (d1.lt.input%structure%epslat) .or. (d2.lt.input%structure%epslat) ) then
             write (*,*)
             write (*, '("Error(fermisurf): zero length plotting vectors")')
             write (*,*)
             stop
           end if
           if ((1.d0-d12 .Lt. input%structure%epslat) .Or. (1.d0+d12 .Lt. input%structure%epslat)) then
              write (*,*)
              write (*, '("Error(fermisurf): zero angle between vectors defining the parallelogram")')
              write (*,*)
              stop
           end if
           nkpt = np2d(1)*np2d(2)
           nkpt_ptr => nkpt
           If (allocated(vkl)) deallocate (vkl)
           Allocate (vkl(3, nkpt))
           vkl_ptr => vkl
           If (allocated(vkc)) deallocate (vkc)
           Allocate (vkc(3, nkpt))
           ik = 0
           Do i1 = 0, np2d(1)-1
              Do i2 = 0, np2d(2)-1
                 ik = ik+1
                 t1 = dble(i1)/dble(np2d(1))
                 t2 = dble(i2)/dble(np2d(2))
                 vkl(:,ik) = t1*vl1(:)+t2*vl2(:)+vclp2d(:,1)
                 Call r3mv(bvec,vkl(:, ik),vkc(:, ik))
              End Do
           End Do

      Else
! determine the k-point grid automatically from radkpt if required
         If (input%groundstate%autokpt) Then
            input%groundstate%ngridk (:) = Int &
           & (input%groundstate%radkpt/&
           & Sqrt(input%structure%crystal%basevect(1, :)**2+&
           & input%structure%crystal%basevect(2, :)**2+&
           & input%structure%crystal%basevect(3, :)**2)) + 1
         End If
! if nktot is set (gt 0), determine the k-point grid automatically from nktot,
! the total number of k-points
         If (input%groundstate%nktot.gt.0) Then
            blen(:)=sqrt(bvec(1,:)**2+bvec(2,:)**2+bvec(3,:)**2)
            lambdab=Dble((input%groundstate%nktot/(blen(1)*blen(2)*blen(3)))**(1.d0/3.d0))
            input%groundstate%ngridk (:) = Max0(1,Int &
           & (lambdab*blen(:)+input%structure%epslat))
             Write (*,*)
             Write (*, '("Info(init1): ngridk determined from nktot: ", 3i8)') &
           & input%groundstate%ngridk(:)
             Write (*,*)
         End If

! setup the default k-point box
         boxl(:, 1) = input%groundstate%vkloff(:) / input%groundstate%ngridk(:)
         boxl(:, 2) = boxl(:, 1)
         boxl(:, 3) = boxl(:, 1)
         boxl(:, 4) = boxl(:, 1)
         boxl(1, 2) = boxl(1, 2) + 1.d0
         boxl(2, 3) = boxl(2, 3) + 1.d0
         boxl(3, 4) = boxl(3, 4) + 1.d0

! allocate the reduced k-point set arrays
         nkptnr = input%groundstate%ngridk(1) * &
         &        input%groundstate%ngridk(2) * &
         &        input%groundstate%ngridk(3)
         If (allocated(ivk)) deallocate (ivk)
         Allocate (ivk(3,nkptnr))
         If (allocated(vkl)) deallocate (vkl)
         Allocate (vkl(3,nkptnr))
         vkl_ptr => vkl
         If (allocated(vkc)) deallocate (vkc)
         Allocate (vkc(3,nkptnr))
         If (allocated(wkpt)) deallocate (wkpt)
         Allocate (wkpt(nkptnr))
         If (allocated(ikmap)) deallocate (ikmap)
         Allocate (ikmap(0:input%groundstate%ngridk(1)-1, &
         &               0:input%groundstate%ngridk(2)-1, &
         &               0:input%groundstate%ngridk(3)-1))

! allocate the non-reduced k-point set arrays
         If (allocated(ivknr)) deallocate (ivknr)
         Allocate (ivknr(3,nkptnr))
         If (allocated(vklnr)) deallocate (vklnr)
         Allocate (vklnr(3,nkptnr))
         If (allocated(vkcnr)) deallocate (vkcnr)
         Allocate (vkcnr(3,nkptnr))
         If (allocated(wkptnr)) deallocate (wkptnr)
         Allocate (wkptnr(nkptnr))
         If (allocated(ikmapnr)) deallocate (ikmapnr)
         Allocate (ikmapnr(0:input%groundstate%ngridk(1)-1, &
         &                 0:input%groundstate%ngridk(2)-1, &
         &                 0:input%groundstate%ngridk(3)-1))

!------------------------------
! generate the k-point set
!------------------------------

         call genppts(.False., .False., &
         &            input%groundstate%ngridk, boxl, nkptnr, &
         &            ikmapnr, ivknr, vklnr, vkcnr, wkptnr)

         call genppts(input%groundstate%reducek, .False., &
         &            input%groundstate%ngridk, boxl, nkpt, &
         &            ikmap, ivk, vkl, vkc, wkpt)
         nkpt_ptr => nkpt

#ifdef TETRA
  ! call to module routine
         If (associated(input%xs)) Then
            If (associated(input%xs%tetra)) Then
               If (input%xs%tetra%tetraocc .Or. tetraopt .Or. input%xs%tetra%tetradf) Then
                  Call genkpts_tet (filext, input%structure%epslat, &
                  &  bvec, maxsymcrys, nsymcrys, lsplsymc, symlat, &
                  &  input%groundstate%reducek, input%groundstate%ngridk, &
                  &  input%groundstate%vkloff, nkpt, ikmap, vkl, wkpt)
               End If
            End If
         End If
#endif

#ifdef XS
  ! determine inverse symmery elements
         Call findsymi(input%structure%epslat, maxsymcrys, nsymcrys, &
         &             symlat, lsplsymc, vtlsymc, isymlat, scimap)
#endif
      End If

      if ( associated(input%groundstate%sirius) ) then
         call setup_sirius_gs_handler(input)
      endif
!
!---------------------!
!     G+k vectors     !
!---------------------!
! determine gkmax
      call Get_gkmax(input, gkmax)

! find the maximum number of G+k-vectors
      if ( associated(input%groundstate%sirius) .and. use_sirius_gkvec ) then
        call get_max_num_gkvec_sirius(ngkmax)
      else
        call getngkmax
      end if

      ngkmax_ptr => ngkmax
! allocate the G+k-vector arrays
      If (allocated(ngk)) deallocate (ngk)
      Allocate (ngk(nspnfv, nkpt))
      ngk_ptr => ngk
      If (allocated(igkig)) deallocate (igkig)
      Allocate (igkig(ngkmax, nspnfv, nkpt))
      If (allocated(vgkl)) deallocate (vgkl)
      Allocate (vgkl(3, ngkmax, nspnfv, nkpt))
      vgkl_ptr => vgkl
      If (allocated(vgkc)) deallocate (vgkc)
      Allocate (vgkc(3, ngkmax, nspnfv, nkpt))
      If (allocated(gkc)) deallocate (gkc)
      Allocate (gkc(ngkmax, nspnfv, nkpt))
      If (allocated(tpgkc)) deallocate (tpgkc)
      Allocate (tpgkc(2, ngkmax, nspnfv, nkpt))
      If (allocated(sfacgk)) deallocate (sfacgk)
      Allocate (sfacgk(ngkmax, natmtot, nspnfv, nkpt))
!      If (allocated(igkfft)) deallocate (igkfft)
!      Allocate (igkfft(ngkmax, nkpt))
!      igkfft=0
      Do ik = 1, nkpt
         Do ispn = 1, nspnfv
            If (isspinspiral()) Then
               sign = (-1.d0)**(ispn+1)
               vl (:) = vkl (:, ik) + sign * 0.5d0 * input%groundstate%spin%vqlss(:)
               vc (:) = vkc (:, ik) + sign * 0.5d0 * vqcss (:)
            Else
               vl (:) = vkl (:, ik)
               vc (:) = vkc (:, ik)
            End If

            if ( associated(input%groundstate%sirius) .and. use_sirius_gkvec ) then
               ! ./seceqn.f90 implies nspnfv=2 is equivalent to spin-spirals
               if ( ispn == 2 ) then
                 write(*,*)"SIRIUS doesn't support spin spirals"
                 call terminate()
               end if
               call get_gkvec_arrays_sirius(ik, ngk(ispn, ik), igkig(:, ispn, ik), vgkl(:, :, ispn, ik), &
                                            vgkc(:, :, ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, ik))
            else
! generate the G+k-vectors
! commented and uncommented versions differ by igkfft(:,ik)
!            Call gengpvec (vl, vc, ngk(ispn, ik), igkig(:, ispn, ik), &
!           & vgkl(:, :, ispn, ik), vgkc(:, :, ispn, ik), gkc(:, ispn, &
!           & ik), tpgkc(:, :, ispn, ik),igkfft(:,ik))
            Call gengpvec (vl, vc, ngk(ispn, ik), igkig(:, ispn, ik), &
           & vgkl(:, :, ispn, ik), vgkc(:, :, ispn, ik), gkc(:, ispn, &
           & ik), tpgkc(:, :, ispn, ik))
           endif
! generate structure factors for G+k-vectors
            Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), &
           & ngkmax, sfacgk(:, :, ispn, ik))
         End Do
      End Do
!
#ifdef XS
      If (init1norealloc) Go To 10
#endif
!---------------------------------!
!     APWs and local-orbitals     !
!---------------------------------!
! allocate linearisation energy arrays
      If (allocated(apwe)) deallocate (apwe)
      Allocate (apwe(maxapword, 0:input%groundstate%lmaxapw, natmtot))
      If (allocated(lorbe)) deallocate (lorbe)
      Allocate (lorbe(maxlorbord, maxlorb, natmtot))
      nlomax = 0
      lolmax = 0
      apwordmax = 0
      Do is = 1, nspecies
! find the maximum APW order
         Do l1 = 0, input%groundstate%lmaxapw
            apwordmax = Max (apwordmax, apword(l1, is))
         End Do
! set the APW linearisation energies to the default
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do l1 = 0, input%groundstate%lmaxapw
               Do io = 1, apword (l1, is)
                  apwe (io, l1, ias) = apwe0 (io, l1, is)
               End Do
            End Do
         End Do
! find the maximum number of local-orbitals
         nlomax = Max (nlomax, nlorb(is))
! set the local-orbital linearisation energies to the default
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ilo = 1, nlorb (is)
               lolmax = Max (lolmax, lorbl(ilo, is))
               Do io = 1, lorbord (ilo, is)
                  lorbe (io, ilo, ias) = lorbe0 (io, ilo, is)
               End Do
            End Do
         End Do
      End Do
      lolmmax = (lolmax+1) ** 2
! generate the local-orbital index
      Call genidxlo
! allocate radial function arrays
      If (allocated(apwfr)) deallocate (apwfr)
      Allocate (apwfr(nrmtmax, 2, apwordmax, &
     & 0:input%groundstate%lmaxapw, natmtot))
      If (allocated(apwdfr)) deallocate (apwdfr)
      Allocate (apwdfr(apwordmax, 0:input%groundstate%lmaxapw, &
     & natmtot))
      If (allocated(lofr)) deallocate (lofr)
      Allocate (lofr(nrmtmax, 2, nlomax, natmtot))
#ifdef XS
10    Continue
#endif
!
!------------------------------------!
!     secular equation variables     !
!------------------------------------!
! number of first-variational states
      nstfv = get_nstfv(chgval, input%groundstate%nempty)
! overlap and Hamiltonian matrix sizes
      If (allocated(nmat)) deallocate (nmat)
      Allocate (nmat(nspnfv, nkpt))
      nmat_ptr => nmat
      If (allocated(npmat)) deallocate (npmat)
      Allocate (npmat(nspnfv, nkpt))
      nmatmax = 0
      nmatmax_ptr => nmatmax
      Do ik = 1, nkpt
         Do ispn = 1, nspnfv
            nmat (ispn, ik) = ngk (ispn, ik) + nlotot
            nmatmax = Max (nmatmax, nmat(ispn, ik))
! packed matrix sizes
            npmat (ispn, ik) = (nmat(ispn, ik)*(nmat(ispn, ik)+1)) / 2
! the number of first-variational states should not exceed the matrix size
            nstfv = Min (nstfv, nmat(ispn, ik))
         End Do
      End Do
! In the standard second-variation implementation the number of basis functions 
! is equal to the first variation state. 
      nbasisfsv = nstfv
! If second variation with local orbitals is used (issvlo=true) the number of &
! second varaition basis is the sum of first variational states and local orbitals. 
      if (issvlo()) then
         nbasisfsv = nstfv + nlotot
      endif
! number of second-variational states
      nstsv = nbasisfsv * nspinor
#ifdef XS
      If (init1norealloc) Go To 20
#endif
! allocate second-variational arrays
      if (allocated(evalsv)) deallocate(evalsv)
      allocate(evalsv(nstsv,nkpt))
      evalsv(:,:) = 0.d0
      if (allocated(occsv)) deallocate(occsv)
      allocate(occsv(nstsv,nkpt))
      occsv(:, :) = 0.d0
      if (allocated(engyknst)) deallocate(engyknst)
      allocate(engyknst(nstfv,nkpt))
      engyknst(:,:) = 0.d0
#ifdef XS
      If (allocated(occsv0)) deallocate (occsv0)
      Allocate (occsv0(nstsv, nkpt))
      occsv (:, :) = 0.d0
      If (allocated(isto0)) deallocate (isto0)
      Allocate (isto0(nkpt))
      isto0 (:) = 0.d0
      If (allocated(isto)) deallocate (isto)
      Allocate (isto(nkpt))
      isto (:) = 0.d0
      If (allocated(istu0)) deallocate (istu0)
      Allocate (istu0(nkpt))
      istu0 (:) = 0.d0
      If (allocated(istu)) deallocate (istu)
      Allocate (istu(nkpt))
      istu (:) = 0.d0
#endif
! allocate overlap and Hamiltonian integral arrays
      If (allocated(oalo)) deallocate (oalo)
      Allocate (oalo(apwordmax, nlomax, natmtot))
      If (allocated(ololo)) deallocate (ololo)
      Allocate (ololo(nlomax, nlomax, natmtot))

        If (allocated(h1aa)) deallocate (h1aa)
        Allocate (h1aa(apwordmax, apwordmax,0:input%groundstate%lmaxapw,natmtot))
        If (allocated(h1loa)) deallocate (h1loa)
        Allocate (h1loa(apwordmax,nlomax,natmtot))
        If (allocated(h1lolo)) deallocate (h1lolo)
        Allocate (h1lolo(nlomax, nlomax, natmtot))
!      endif
! allocate and generate complex Gaunt coefficient array
      If (allocated(gntyry)) deallocate (gntyry)
      Allocate (gntyry(lmmaxapw, lmmaxapw, lmmaxapw))
      nonzcount=0
      Do l1 = 0, input%groundstate%lmaxapw
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do l2 = 0, input%groundstate%lmaxapw
               Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        gntyry (lm1, lm2, lm3) = gauntyry (l1, l2, l3, m1, m2, m3)
                        if ((abs(gntyry (lm1, lm2, lm3)).gt.1d-20)) nonzcount=nonzcount+1
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      If (allocated(gntryy)) deallocate (gntryy)
      If (allocated(gntnonz)) deallocate (gntnonz)
      If (allocated(gntnonzlm1)) deallocate (gntnonzlm1)
      If (allocated(gntnonzlm2)) deallocate (gntnonzlm2)
      If (allocated(gntnonzlm3)) deallocate (gntnonzlm3)
      If (allocated(gntnonzlindex)) deallocate (gntnonzlindex)
      If (allocated(gntnonzl2index)) deallocate (gntnonzl2index)
      Allocate (gntryy(lmmaxvr, lmmaxapw, lmmaxapw))
      allocate(gntnonz(nonzcount),gntnonzlm1(nonzcount+1),gntnonzlm2(nonzcount),gntnonzlm3(nonzcount+1))
      allocate(gntnonzlindex(0:input%groundstate%lmaxmat))
      allocate(gntnonzl2index(lmmaxmat,lmmaxmat))

      Do l1 = 0, input%groundstate%lmaxapw
        Do m1 = - l1, l1
          lm1 = idxlm (l1, m1)
          Do l3 = 0, input%groundstate%lmaxapw
            Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)
              Do l2 = 0, input%groundstate%lmaxvr
                Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  gntryy (lm2, lm3, lm1) = gauntyry (l1, l2, l3, m1, m2, m3)
                End Do
              End Do
            End Do
          End Do
        End Do
      End Do

      i1=0
      Do l1 = 0, input%groundstate%lmaxmat
        gntnonzlindex(l1)=i1+1
        Do m1 = - l1, l1
          lm1 = idxlm (l1, m1)
          Do l3 = 0, input%groundstate%lmaxmat
            Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)
              gntnonzl2index(lm1,lm3)=i1+1
              Do l2 = 0, input%groundstate%lmaxvr
                Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  if ((abs(gntryy (lm2, lm3, lm1)).gt.1d-20).and.(lm1.le.lmmaxmat).and.(lm3.le.lmmaxmat)) then
                    i1=i1+1
                    gntnonz(i1)=gntryy(lm2, lm3, lm1)
                    gntnonzlm1(i1)=lm1
                    gntnonzlm3(i1)=lm3
                    gntnonzlm2(i1)=lm2
                  endif
                End Do
              End Do
            End Do
          End Do
        End Do
      End Do
      gntnonzlm3(nonzcount+1)=0
      gntnonzlm1(nonzcount+1)=0

      ! compact Gaunt coeffiecient array
      if (allocated(indgnt)) deallocate(indgnt)
      if (allocated(listgnt)) deallocate(listgnt)
      allocate(indgnt(input%groundstate%lmaxvr+1,lmmaxapw,lmmaxapw))
      allocate(listgnt(input%groundstate%lmaxvr+1,lmmaxapw,lmmaxapw))
      indgnt(:,:,:)=0
      listgnt(:,:,:)=zzero
      maxi=0
      do lm3=1,lmmaxapw
        do lm1=1,lmmaxapw
          i=0
          do lm2=1,(input%groundstate%lmaxvr+1)**2 !lmmaxvr
            if ((dble(gntryy(lm2,lm1,lm3)).ne.0d0).or.(dimag(gntryy(lm2,lm1,lm3)).ne.0d0)) then
              i=i+1
              listgnt(i,lm1,lm3)=gntryy(lm2,lm1,lm3)
              indgnt(i,lm1,lm3)=lm2
            endif
          enddo
          maxi=max(i,maxi)
        enddo
      enddo

#ifdef XS
20    Continue
#endif
!
      nullify(arpackinverse)
      Call timesec (ts1)
      call stopwatch("exciting:init1", 0)
!
      Return
End Subroutine
!EOC
