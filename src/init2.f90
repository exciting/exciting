
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine init2
      Use modmain
      Use modinput
#ifdef XS
      Use modxs
#endif
      Implicit None
! local variables
      logical :: redq
      Integer :: is, ia, ist, ic, m
      Real (8) :: ts0, ts1
      Real (8) :: boxl (3, 4)
#ifdef XS
      Real (8) :: v (3), t1
      Integer :: iq, iv (3)
      Character (256) :: filex
#endif

      Call timesec (ts0)

!---------------------!
!     q-point set     !
!---------------------!
      redq = .true.
      ! phonons
      if ((task .eq. 200).or.(task .eq. 201).or. &
          (task .eq. 210).or.(task .eq. 220).or.(task .eq. 230).or. &
          (task .eq. 240).or.(task .eq. 245).or.(task .eq. 250)) then
          ngridq(:)=input%phonons%ngridq(:)
          redq=input%phonons%reduceq
      end if
      ! phonons for q-points from list
      if (task .eq. 230) then
          nphwrt = size (input%phonons%qpointset%qpoint, 2)
      	  if (allocated(vqlwrt)) deallocate(vqlwrt)
          allocate(vqlwrt(3,nphwrt))
          Do iq = 1, nphwrt
            vqlwrt(:,iq) = input%phonons%qpointset%qpoint(:, iq)
          end do
      end if
      ! OEP, Hartree-Fock or RDMFT
      If (associated(input%groundstate%OEP) .or. &
      &   associated(input%groundstate%Hybrid) .or. &
      &  (task .Eq. 300) .or. &
      &  (associated(input%groundstate%HartreeFock))) Then 
        ngridq (:) = input%groundstate%ngridk(:)
        redq = .False.
      End If

#ifdef XS
      If (task .Le. 300) Then
#endif
! allocate the q-point arrays
         If (allocated(ivq)) deallocate (ivq)
         Allocate (ivq(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vql)) deallocate (vql)
         Allocate (vql(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vqc)) deallocate (vqc)
         Allocate (vqc(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(wqpt)) deallocate (wqpt)
         Allocate (wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(iqmap)) deallocate (iqmap)
         Allocate (iqmap(0:ngridq(1)-1, 0:ngridq(2)-1, 0:ngridq(3)-1))
! setup the q-point box (offset should always be zero)
         boxl (:, :) = 0.d0
         boxl (1, 2) = 1.d0
         boxl (2, 3) = 1.d0
         boxl (3, 4) = 1.d0
! generate the q-point set, note that the vectors vql and vqc are mapped to the
! first Brillouin zone
         Call genppts(redq, .True., ngridq, boxl, &
         &            nqpt, iqmap, ivq, vql, vqc, wqpt)
#ifdef XS
      End If
#endif
#ifdef XS

! Q-/q-point set should have no offset
! setup the q-point box (offset should always be zero)
      boxl(:,:) = 0.d0
      boxl(1,2) = 1.d0
      boxl(2,3) = 1.d0
      boxl(3,4) = 1.d0
      ! boxl(:,2) = boxl(:,2) + boxl(:,1)
      ! boxl(:,3) = boxl(:,3) + boxl(:,1)
      ! boxl(:,4) = boxl(:,4) + boxl(:,1)

! assign momentum transfer Q-points set to q-point set
      If ((task .Ge. 301) .And. (task .Le. 399)) Then
         nqpt = size (input%xs%qpointset%qpoint, 2)
         If (allocated(vqlmt)) deallocate (vqlmt)
         Allocate (vqlmt(3, nqpt))
         If (allocated(ivgmt)) deallocate (ivgmt)
         Allocate (ivgmt(3, nqpt))
         If (allocated(vql)) deallocate (vql)
         Allocate (vql(3, nqpt))
         If (allocated(vqc)) deallocate (vqc)
         Allocate (vqc(3, nqpt))
         Do iq = 1, nqpt
            v(:) = input%xs%qpointset%qpoint(:, iq)
            iv(:) = 0
            ! map Q-point to reciprocal unit cell
            If (input%xs%tddft%mdfqtype .Eq. 1) Call r3frac(input%structure%epslat, v, iv)
            vqlmt(:,iq) = v(:)
            ivgmt(:,iq) = iv(:)
            vql(:,iq) = vqlmt(:,iq)
            vqc(:,iq) = vql(1,iq)*bvec(:,1) + &
            &           vql(2,iq)*bvec(:,2) + &
            &           vql(3,iq)*bvec(:,3)
            ! check consistency of Q-point with gqmax
            v(:) = input%xs%qpointset%qpoint(1,iq)*bvec(:,1)+ &
            &      input%xs%qpointset%qpoint(2,iq)*bvec(:,2)+ &
            &      input%xs%qpointset%qpoint(3,iq)*bvec(:,3)
            t1 = sqrt(v(1)**2+v(2)**2+v(3)**2)
            if ((input%xs%gqmax.ne.0.d0).and.(t1.gt.input%xs%gqmax)) then
              write(*,*)
              write(*,'("Info(init2): Q-point exceeds gqmax")')
              write(*,'(" Q-point number         : ",i6)') iq
              write(*,'(" Q-point (latt. coords.): ",3g18.10)') vql(:,iq)
              write(*,'(" Q-point length         : ",3g18.10)') t1
              write(*,'(" gqmax                  : ",g18.10)') input%xs%gqmax
              write(*,*)
              if (input%xs%gqmaxtype .eq. "|G+q|") then
                write(*,'("Error(init2): gqmaxtype = ""|G+q|"" and Q-point exceeds gqmax")')
                write(*,'(" set gqmaxtype to ""|G|""")')
                write(*,*)
                stop
              end if
            end if
            if (input%xs%tddft%mdfqtype .eq. 1) then
              write(*,'("Error(init2): mdfqtype = 1 and Q-point exceeds gqmax")')
              write(*,'(" set mdfqtype to 0")')
              write(*,*)
              stop
            end if
            ! check if required tasks can be managed for non-zero Q-point
            if (t1.ne.0.d0) then
              if (input%xs%xstype.eq."TDDFT") then
                if ((input%xs%xstype.eq."MB1").or.(input%xs%xstype.eq."MB1_NLF")) then
                  write(*,*)
                  write(*,'("Error(init2): BSE derived xc kernel only works for optics (Q=0) - code limitation")')
                  write(*,*)
                  stop
                end if
              end if
              if (input%xs%xstype.eq."BSE") then
                write(*,*)
                write(*,'("Error(init2): BSE only works for optics (Q=0) - code limitation")')
                write(*,*)
                stop
              end if
            end if
         End Do
      Else if (task .ge. 400) then
         ! determine only integer-part of Q-points
         If (allocated(ivgmt)) deallocate (ivgmt)
         Allocate (ivgmt(3, size(input%xs%qpointset%qpoint, 2)))
         Do iq = 1, size(input%xs%qpointset%qpoint, 2)
            v(:) = input%xs%qpointset%qpoint(:,iq)
            iv(:) = 0
            ! map Q-point to reciprocal unit cell
            If (input%xs%tddft%mdfqtype .Eq. 1) Call r3frac(input%structure%epslat, v, iv)
            ivgmt(:,iq) = iv(:)
         End Do
      End If

      ! generate q-point set from grid
      If ((task .Ge. 400) .And. (task .Le. 439)) Then
         If (allocated(ivq)) deallocate (ivq)
         Allocate (ivq(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vql)) deallocate (vql)
         Allocate (vql(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vqc)) deallocate (vqc)
         Allocate (vqc(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(wqpt)) deallocate (wqpt)
         Allocate (wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(iqmap)) deallocate (iqmap)
         Allocate (iqmap(0:ngridq(1)-1, 0:ngridq(2)-1, 0:ngridq(3)-1))
         ! generate reduced q-point set
         Call genppts(input%xs%reduceq, input%xs%BSE%fbzq, ngridq, &
         &            boxl, nqpt, iqmap, ivq, vql, vqc, wqpt)
         nqptr = nqpt
      End If
      If ((task .Eq. 440) .Or. (task .Eq. 441) .Or. &
      &   (task .Eq. 445) .Or. (task .Eq. 446) .Or. &
      &   (task .Eq. 450) .Or. (task .Eq. 451) .Or. &
      &   (task .Eq. 499) .Or. (task .Eq. 700) .Or. &
      &   (task .Eq. 710)) Then
         If (allocated(ivqr)) deallocate (ivqr)
         Allocate (ivqr(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vqlr)) deallocate (vqlr)
         Allocate (vqlr(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vqcr)) deallocate (vqcr)
         Allocate (vqcr(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(wqptr)) deallocate (wqptr)
         Allocate (wqptr(ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(iqmapr)) deallocate (iqmapr)
         Allocate (iqmapr(0:ngridq(1)-1, 0:ngridq(2)-1, 0:ngridq(3)-1))
         ! generate reduced q-point set
         Call genppts(input%xs%reduceq, input%xs%BSE%fbzq, ngridq, &
         &            boxl, nqptr, iqmapr, ivqr, vqlr, vqcr, wqptr)
         If (allocated(ivq)) deallocate (ivq)
         Allocate (ivq(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vql)) deallocate (vql)
         Allocate (vql(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(vqc)) deallocate (vqc)
         Allocate (vqc(3, ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(wqpt)) deallocate (wqpt)
         Allocate (wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
         If (allocated(iqmap)) deallocate (iqmap)
         Allocate (iqmap(0:ngridq(1)-1, 0:ngridq(2)-1, 0:ngridq(3)-1))
        ! generate non-reduced q-point set
         Call genppts(.False., input%xs%BSE%fbzq, ngridq, boxl, nqpt, &
         &            iqmap, ivq, vql, vqc, wqpt)
      End If
! find (little/small) group of q
      If (allocated(nsymcrysq)) deallocate (nsymcrysq)
      Allocate (nsymcrysq(nqpt))
      If (allocated(scqmap)) deallocate (scqmap)
      Allocate (scqmap(nsymcrys, nqpt))
      If (allocated(ivscwrapq)) deallocate (ivscwrapq)
      Allocate (ivscwrapq(3, nsymcrys, nqpt))
!do iq=1,nqpt
!   call findgroupq(fbzq,vql(1,iq),epslat,bvec,symlat,nsymcrys,lsplsymc,&
!	nsymcrysq(iq),scqmap(1,iq),ivscwrapq(1,1,iq))
!end do
  if (task .ge. 301) then

!-----------------------!
!     k+q-point set	    !
!-----------------------!
      If (allocated(qvkloff)) deallocate (qvkloff)
      Allocate (qvkloff(3, 0:nqpt))
      If (allocated(ikmapikq)) deallocate (ikmapikq)
      Allocate (ikmapikq(nkpt, nqpt))
      qvkloff (:, 0) = input%groundstate%vkloff(:)
      Do iq = 1, nqpt
        ! offset for k+q-point set derived from q-point
        Call genqvkloff(vql(1,iq), qvkloff(1,iq))
        ! map from k-point index to k+q point index for same k
        Call findkmapkq(vql(1,iq), qvkloff(1,iq), ikmapikq(1,iq))
      End Do
!
!---------------------!
!     G+q-point set   !
!---------------------!
! warning for small gqmax
      If (input%xs%gqmax .Ge. gkmax) Then
        write(*,'(a, 2g18.10)') 'Warning(init2/xs): input%xs%gqmax >= gkmax: ', &
        &  input%xs%gqmax, gkmax
      End If
! check consistency with FFT G-vector array
      if (input%groundstate%gmaxvr .lt. 2*gkmax + input%xs%gqmax) then
         write(*,*)
         write(*,'("Error(init2): gmaxvr < 2 gkmax + gqmax : ",2g18.10)') &
           & input%groundstate%gmaxvr, 2*gkmax+input%xs%gqmax
         write(*,'(" gmaxvr : ",g18.10)') input%groundstate%gmaxvr
         write(*,'(" gkmax  : ",g18.10)') gkmax
         write(*,'(" gqmax  : ",g18.10)') input%xs%gqmax
         write(*,'(" maximum value for gqmax    : ",g18.10)') &
           & input%groundstate%gmaxvr - 2*gkmax
         write(*,'(" Increase gmaxvr and re-do SCF calculation or decrease gqmax or rgkmax.")')
         write(*,*)
         call terminate
      end if
! maximum number of G+q vectors for all q
      Call getngqmax
! allocate the G+q-vector arrays
      If (allocated(ngq)) deallocate (ngq)
      Allocate (ngq(nqpt))
      If (allocated(igqig)) deallocate (igqig)
      Allocate (igqig(ngqmax, nqpt))
      If (allocated(vgql)) deallocate (vgql)
      Allocate (vgql(3, ngqmax, nqpt))
      If (allocated(vgqc)) deallocate (vgqc)
      Allocate (vgqc(3, ngqmax, nqpt))
      If (allocated(gqc)) deallocate (gqc)
      Allocate (gqc(ngqmax, nqpt))
      If (allocated(tpgqc)) deallocate (tpgqc)
      Allocate (tpgqc(2, ngqmax, nqpt))
      If (allocated(sfacgq)) deallocate (sfacgq)
      Allocate (sfacgq(ngqmax, natmtot, nqpt))
      If (allocated(ylmgq)) deallocate (ylmgq)
      Allocate (ylmgq(lmmaxapw, ngqmax, nqpt))
      If (allocated(ivgigq)) deallocate (ivgigq)
      Allocate (ivgigq(intgqv(1,1):intgqv(1,2), &
      &                intgqv(2,1):intgqv(2,2), &
      &                intgqv(3,1):intgqv(3,2), nqpt))
      Do iq = 1, nqpt
   ! generate G+q vectors
         Call gengqvec(iq, vql(1,iq), vqc(1,iq), ngq(iq), &
         &             igqig(1,iq), vgql(1,1,iq), vgqc(1,1,iq), &
         &             gqc(1,iq), tpgqc(1,1,iq))
   ! generate structure factors for G-vectors
         Call gensfacgp(ngq(iq), vgqc(1,1,iq), ngqmax, sfacgq(1,1,iq))
   ! spherical harmonics for G+q-vectors
         Call genylmgq(iq, input%groundstate%lmaxvr)
      End Do
!
!---------------------------!
!     Coulomb potential     !
!---------------------------!
      If (allocated(sptclg)) deallocate (sptclg)
      Allocate (sptclg(ngqmax, nqpt))
      Do iq = 1, nqpt
   ! set up Coulomb potential square root
         Call genptclg ('nocutoff', ngqmax, ngq(iq), vgqc(:, :, iq), &
        & gqc(:, iq), sptclg(:, iq))
      End Do
!
!------------------------!
!     radial functions   !
!------------------------!
! read density and potentials from file (STATE.OUT) exclusively
      isreadstate0 = .True.
      If (hybridhf) Then
! in case of HF hybrids read PBE potential
            isreadstate0 = .false. 
            filex=filext
            filext='_PBE.OUT'
            Call readstate
            filext=filex
      Else If (input%xs%dogroundstate .Ne. "fromscratch") Then 
         Call readstate
      Else
         If(task .Ne. 301) Then
            isreadstate0 = .False.
            filex = trim(filext)
            filext = '_QMT001.OUT'
            Call readstate
            filext = trim(filex)
         End If
      End If
      isreadstate0 = .False.
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! update potential in case of HF Hybrids
        If (hybridhf)  Then
            filex= trim(filext)
            filext='.OUT'
            Call readstate
            filext=trim(filex)
        End If

  ! end for task >= 301 case
  end if
#endif
!
!-----------------------------------------------!
!     OEP, Hartree-Fock and RDMFT variables     !
!-----------------------------------------------!
      If ((input%groundstate%xctypenumber .Lt. 0) .Or. (task .Eq. 5) &
     & .Or. (task .Eq. 6) .Or. (task .Eq. 300) .Or.  (xctype(2) .Ge. 400)&
     &.Or.  (xctype(1) .Ge. 400)) Then
! determine the 1/q^2 integral weights if required
         Call genwiq2
! output the 1/q^2 integrals to WIQ2.OUT
         Call writewiq2
      End If
!      If ((input%groundstate%xctypenumber .Lt. 0) .Or.  (xctype(2) .Ge. 400)&
!     &.Or.  (xctype(1) .Ge. 400)) Then
      If (associated(input%groundstate%OEP)) Then
! initialise OEP residual magnitude
         resoep = 1.d0
! find maximum core states over all species
         ncrmax = 0
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ic = 0
               Do ist = 1, spnst (is)
                  If (spcore(ist, is)) Then
                     Do m = - spk (ist, is), spk (ist, is) - 1
                        ic = ic + 1
                     End Do
                  End If
               End Do
               ncrmax = Max (ncrmax, ic)
            End Do
         End Do
! allocate and zero the complex exchange potential and field
         If (allocated(zvxmt)) deallocate (zvxmt)
         Allocate (zvxmt(lmmaxvr, nrcmtmax, natmtot))
         zvxmt (:, :, :) = 0.d0
         If (allocated(zvxir)) deallocate (zvxir)
         Allocate (zvxir(ngrtot))
         zvxir (:) = 0.d0
         If (associated(input%groundstate%spin)) Then
            If (allocated(zbxmt)) deallocate (zbxmt)
            Allocate (zbxmt(lmmaxvr, nrcmtmax, natmtot, ndmag))
            zbxmt (:, :, :, :) = 0.d0
            If (allocated(zbxir)) deallocate (zbxir)
            Allocate (zbxir(ngrtot, ndmag))
            zbxir (:, :) = 0.d0
         End If
      End If
      If ((task .Eq. 5) .Or. (task .Eq. 6) .Or. (task .Eq. 300)) Then
! allocate the kinetic matrix elements for Hartree-Fock/RDMFT
         If (allocated(kinmatc)) deallocate (kinmatc)
         Allocate (kinmatc(nstsv, nstsv, nkpt))
      End If
      If (task .Eq. 300) Then
         If (allocated(vclmat)) deallocate (vclmat)
         Allocate (vclmat(nstsv, nstsv, nkpt))
         If (allocated(dkdc)) deallocate (dkdc)
         Allocate (dkdc(nstsv, nstsv, nkpt))
         If (allocated(vnlrdm)) deallocate (vnlrdm)
         Allocate (vnlrdm(nstsv, nkpt, nstsv, nkptnr))
      End If
!
      Call timesec (ts1)
!      timeinit = timeinit + ts1 - ts0
!
      Return
End Subroutine
