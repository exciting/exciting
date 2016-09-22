! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bandstr
! !INTERFACE:
!
!
Subroutine bandstr
  ! !USES:
  Use modinput
  Use modmain
  Use modmpi
  Use FoX_wxml

  ! !DESCRIPTION:
  !   Produces a band structure along the path in reciprocal-space which connects
  !   the vertices in the array {\tt vvlp1d}. The band structure is obtained from
  !   the second-variational eigenvalues and is written to the file {\tt BAND.OUT}
  !   with the Fermi energy set to zero. If required, band structures are plotted
  !   to files {\tt BAND\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species {\tt ss},
  !   which include the band characters for each $l$ component of that atom in
  !   columns 4 onwards. Column 3 contains the sum over $l$ of the characters.
  !   Vertex location lines are written to {\tt BANDLINES.OUT}.
  !
  ! !REVISION HISTORY:
  !   Created June 2003 (JKD)
  !   Modified June 2012 (DIN)
  !   Modified March 2014 (UW)
  !EOP
  !BOC
  Implicit None
  ! local variables
  Integer :: lmax, lmmax, l, m, lm
  Integer :: ik, ispn, is, ia, ias, iv, ist
  Real (8) :: emin, emax, sum
  Character (256) :: fname
  ! allocatable arrays
  Real (8), Allocatable :: evalfv (:, :)
  ! low precision for band character array saves memory
  Real (4), Allocatable :: bc (:, :, :, :)
  Complex (8), Allocatable :: dmat (:, :, :, :, :)
  Complex (8), Allocatable :: apwalm (:, :, :, :, :)
  Complex (8), Allocatable :: evecfv (:, :, :)
  Complex (8), Allocatable :: evecsv (:, :)
  Character (128) :: buffer, buffer1
  Type (xmlf_t), Save :: xf
  !
  integer :: nkpt0, nstsv0, Recl
  real(8),    Allocatable :: vkl0(:,:), ehf(:,:), deltax(:,:)
  complex(8), allocatable :: e0(:,:), e1(:,:)
  logical :: exist,hybcheck
!WANNIER
  integer :: ir, iq, ig, io, ilo, iproj, iband
  integer :: nr, nproj, nrpt, niband, npband, pbandstart, pbandstop
  integer :: isproj, iaproj, iasproj, iloproj
  real (8) :: r2, u, f
  real (8), allocatable :: uf(:), gr(:), cf(:,:)
  real (8), allocatable :: rolpi(:,:), oaa(:,:,:,:), rpt(:,:)
  integer, allocatable :: projst(:,:), idxnn(:,:)
  complex (8), allocatable :: projme(:,:), olpme(:,:), transme(:,:,:), wanme(:,:,:), auxmat(:,:), auxmat2(:,:)
  complex (8) :: auxc
  real (8), allocatable :: evaltmp(:), vklcpy(:,:), vkccpy(:,:)
  complex (8), allocatable :: mattmp(:,:), evectmp(:,:), evalcplx(:)

  complex (8), allocatable :: lapack_work(:), lapack_rwork(:)
  integer :: lapack_lwork, lapack_info
  integer, allocatable :: lapack_iwork(:)

  complex (8), allocatable :: mzero(:,:,:,:)
  integer :: ngridk(3)

  complex(8), allocatable :: mlwf0(:,:,:,:)
!END WANNIER
! initialise universal variables
  Call init0
  Call init1

if (input%properties%bandstructure%wannier) then
  !--------------------------------------------------!      
  ! Calculate bandstructure by Wannier interpolation !
  !--------------------------------------------------!
  nrpt = nkpt

  ! maximum angular momentum for band character
  lmax = min( 3, input%groundstate%lmaxapw)
  lmmax = (lmax+1)**2
  if (input%properties%bandstructure%character) then
     allocate( bc( 0:lmax, natmtot, nstsv, nkpt))
  end If
  ! read density and potentials from file
  Call readstate
  ! read Fermi energy from file
  Call readfermi                !saves fermi energy in variable 'efermi'
  ! find the new linearisation energies
  Call linengy
  ! generate the APW radial functions
  Call genapwfr                 !generates APW radial functions and save them in 'apwfr[5]'
  ! generate the local-orbital radial functions
  Call genlofr                  !generates LO radial functions and save them in 'lofr[4]'
  Call olprad

  allocate( apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  allocate( evalfv(nstfv, nspnfv), evecfv(nmatmax, nstfv, nspnfv))
  allocate( uf(nrmtmax), gr(nrmtmax), cf(3, nrmtmax), idxnn( nkpt, 6))
  
  !============================================
  ! read projection local-orbitals
  allocate( projst(nlotot, 5))      ! (index of projector ; 1=index of species, 2=index of atom, 3=index of LO within species, 4=angular momentum, 5=magnetic quantum number)
  projst(:,:) = 0
  nproj = 0
  DO is = 1, nspecies
    DO ia = 1, natoms( is)
      DO ilo = 1, nlorb( is)
        IF( lorbwfproj( ilo, is)) THEN
          l = lorbl( ilo, is)
          DO m = -l, l
            nproj = nproj + 1
            projst( nproj, 1) = is
            projst( nproj, 2) = ia
            projst( nproj, 3) = ilo
            projst( nproj, 4) = l
            projst( nproj, 5) = m 
          END DO
        END IF
      END DO
    END DO
  END DO

  ! write next neighbours to array
  ngridk = input%groundstate%ngridk
  DO ir = 1, ngridk( 3)
    DO ia = 1, ngridk( 2)
      DO io = 1, ngridk( 1)
        ik = modulo( ir-1, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia-1, ngridk( 2))*ngridk( 1) + modulo( io-1, ngridk(1))+1
        idxnn( ik, 1) = modulo( ir-1, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia-1, ngridk( 2))*ngridk( 1) + modulo( io, ngridk(1))+1
        idxnn( ik, 2) = modulo( ir-1, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia, ngridk( 2))*ngridk( 1) + modulo( io-1, ngridk(1))+1
        idxnn( ik, 3) = modulo( ir, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia-1, ngridk( 2))*ngridk( 1) + modulo( io-1, ngridk(1))+1
        idxnn( ik, 4) = modulo( ir-1, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia-1, ngridk( 2))*ngridk( 1) + modulo( io-2, ngridk(1))+1
        idxnn( ik, 5) = modulo( ir-1, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia-2, ngridk( 2))*ngridk( 1) + modulo( io-1, ngridk(1))+1
        idxnn( ik, 6) = modulo( ir-2, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
            modulo( ia-1, ngridk( 2))*ngridk( 1) + modulo( io-1, ngridk(1))+1
      END DO
    END DO
  END DO
                                              !
  !--------------------------------------------

  !============================================
  ! allocating some variables
  niband = 4 
  pbandstart = 1 
  pbandstop = 4
  npband = pbandstop - pbandstart + 1
  niband = MIN( niband, nproj)
  nproj = MIN( niband, nproj)
  allocate( rolpi( apwordmax, nproj))                   ! radial overlap integrals
  allocate( projme( npband, nproj))                     ! projection matrix elements	-- A(i,j)
  allocate( olpme( nproj, nproj))                       ! overlap matrix 		-- S(i,j) 	/ S = adj(A)*A
  allocate( transme( nkpt, npband, nproj))              ! transformation matrix		-- V(k,i,j)	/ V = A*S^(-1/2)
  allocate( wanme( nrpt, nproj, nproj))                 ! Wannier matrix elements       -- <0i|H|Rj>
  allocate( evaltmp( nproj), evectmp( nproj, nproj))    ! temporary eigenvalues and eigenvectors for matrix diagonalization
  allocate( auxmat( nproj, nproj))                      ! maxtrix container for auxiliary matrices for calculational purposes
  allocate( rpt( 3, nrpt))                              ! set of R-points

  allocate( mzero( npband, npband, nkpt, 6))
                                              !
  !--------------------------------------------
  write(*,'(A)', advance='no') "bands "
  write(*,'(I2)', advance='no') pbandstart
  write(*,'(A)', advance='no')  " to " 
  write(*,'(I2)', advance='no') pbandstop
  write(*,'(A)', advance='no') " are used for projection ("
  write(*,'(I2)', advance='no') npband
  write(*,*) " bands)"
  write(*,*) nproj, " projection functions"
  write(*,*) niband, " bands are going to be interpolated"
 
  !============================================ 
  ! radial overlap integrals
  rolpi(:,:) = zzero
  DO iproj = 1, nproj
    uf(:) = 0d0
    is = projst( iproj, 1)
    nr = nrmt( is)
    ias = idxas( projst( iproj, 2), is)
    lm = idxlm( projst( iproj, 4), projst( iproj, 5))
    DO io = 1, apword( projst( iproj, 4), is)
      DO ir = 1, nr
        r2 = spr( ir, is)**2
        u = apwfr( ir, 1, io, projst( iproj, 4), ias)
        f = lofr( ir, 1, projst( iproj, 3), ias)
        uf( ir) = u*f*r2
      END DO
      Call fderiv(-1, nr, spr(:, is), uf, gr, cf)
      rolpi( io, iproj) = gr( nr)
    END DO
  END DO
                                              !
  !--------------------------------------------

  DO ik = 1, nkpt
    !============================================
    ! get eigenenergies, basis function coefficients and matching coefficients
    Call getevecfv( vkl(:, ik), vgkl(:, :, :, ik), evecfv)
    Call getevalfv( vkl(:, ik), evalfv)
    Call match( ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, &
         & 1, ik), sfacgk(:, :, 1, ik), apwalm(:, :, :, :, &
         & 1))
                                                !
    !--------------------------------------------
    !write(*,*) mod( ik-1, 8), mod(( ik-1)/8, 8), (ik-1)/64, evalfv( 4, 1) 
    !write(*,'(F12.10,1x)', advance='no') evalfv( 4, 1)-efermi

    !============================================  
    ! projection matrix elements and overlap matrix elements
    DO iproj = 1, nproj
      is = projst( iproj, 1)
      nr = nrmt( is)
      ias = idxas( projst( iproj, 2), is)
      lm = idxlm( projst( iproj, 4), projst( iproj, 5))
      DO iband = 1, npband
        projme( iband, iproj) = zzero
        DO ig = 1, ngk( 1, ik)
          auxc = zzero
          DO io = 1, apword( projst( iproj, 4), is)
            auxc = auxc + DCONJG( apwalm( ig, io, lm, ias, 1))*rolpi( io, iproj)
          END DO
          projme( iband, iproj) = projme( iband, iproj) + DCONJG( evecfv( ig, pbandstart+iband-1, 1))*auxc
        END DO
        DO io = 1, nlorb( is)
          l = lorbl( io, is)
          IF( l .EQ. projst( iproj, 4)) THEN
            projme( iband, iproj) = projme( iband, iproj) + DCONJG( evecfv( ngk( 1, ik)+idxlo( lm, io, ias), pbandstart+iband-1, 1))*ololo( io, projst( iproj, 3), ias)
          END IF
        END DO
      END DO
    END DO
    olpme(:,:) = zzero
    Call ZGEMM( 'C', 'N', nproj, nproj, npband, (1d0, 0d0), projme, npband, projme, npband, zzero, olpme, nproj)
                                                !
    !--------------------------------------------

    !============================================    
    ! transformation matrix
    Call diaghermat( nproj, olpme, evaltmp, evectmp)
    auxmat(:,:) = zzero
    DO io = 1, nproj
      IF( evaltmp( io).EQ.0d0) THEN
        write(*,*) "0 eigenvalue ", ik, io
      ELSE
        auxmat( io, io) = 1d0/sqrt( evaltmp( io)) 
      END IF
    END DO

    Call ZGEMM( 'N', 'N', nproj, nproj, nproj, (1d0, 0d0), evectmp, nproj, auxmat, nproj, zzero, olpme, nproj)
    Call ZGEMM( 'N', 'C', nproj, nproj, nproj, (1d0, 0d0), olpme, nproj, evectmp, nproj, zzero, auxmat, nproj)
    Call ZGEMM( 'N', 'N', npband, nproj, nproj, (1d0, 0d0), projme, npband, auxmat, nproj, zzero, transme( ik, :, :), npband)
                                                !
    !--------------------------------------------
  END DO
  Deallocate( projme, olpme, auxmat, evaltmp, evectmp)
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  allocate( mlwf0( npband, npband, 6, nkpt))
  do ik = 1, nkpt
    do iq = 1, 6
      Call emat_wannier( ik, vkl( :, idxnn( ik, iq)), pbandstart, npband, &
      & pbandstart, npband, mlwf0( :, :, iq, ik)) 
      write(*,*) mlwf0( :, :, iq, ik)
    end do
  end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !============================================
    ! Wannier matrix

  ! read lattice site grid
  ia = 0
  DO ir = -input%groundstate%ngridk(3)/2, -input%groundstate%ngridk(3)/2+input%groundstate%ngridk(3)-1
    DO io = -input%groundstate%ngridk(2)/2, -input%groundstate%ngridk(2)/2+input%groundstate%ngridk(2)-1
      DO is = -input%groundstate%ngridk(1)/2, -input%groundstate%ngridk(1)/2+input%groundstate%ngridk(1)-1
        ia = ia + 1
        rpt( :, ia) = ir*input%structure%crystal%basevect( :, 3) + io*input%structure%crystal%basevect( :, 2) + is*input%structure%crystal%basevect( :, 1) 
      END DO
    END DO
  END DO

  ! calculate Hamlitonian matrix elements in Wannier representation 
  Allocate( auxmat( nproj, nproj))
  DO ir = 1, nrpt 
    wanme( ir, :, :) = zzero
    DO ia = 1, nproj
      DO is = 1, nproj
        DO ik = 1, nkpt
          auxc = zzero
          Call getevalfv( vkl( :, ik), evalfv)
          DO iband = 1, npband 
            auxc = auxc + evalfv( pbandstart+iband-1, 1)*transme( ik, iband, is)*DCONJG( transme( ik, iband, ia))
          END DO
          wanme( ir, ia, is) = wanme( ir, ia, is) + exp( -zi*DOT_PRODUCT( rpt( :, ir), vkc( :, ik)))*auxc
        END DO
      END DO
    END DO
    wanme( ir, :, :) = wanme( ir, :, :)/nkpt
    !write(*,*) wanme( ir, :, :)
  END DO
  Deallocate( transme, auxmat)
  
                                                ! 
    !--------------------------------------------

  IF( allocated( vklcpy)) Deallocate( vklcpy)
  Allocate( vklcpy( 3, nkpt))
  vklcpy = vkl
  IF( allocated( cf)) Deallocate( cf)
  Allocate( cf( nkpt, nstfv))
  DO ik = 1, nkpt
    Call getevalfv( vkl( :, ik), evalfv)
    cf( ik, :) = evalfv( :, 1)
  END DO
    
    !============================================
    ! interpolation
  write(*,*) "start interpolation"
  input%properties%bandstructure%wannier = .FALSE. 
  Call init1        ! initialize k-points along path
  Allocate( auxmat( nproj, nproj), evectmp( nproj, nproj))
  if( allocated( evalsv)) Deallocate( evalsv)
  Allocate( evalsv( nproj, nkpt), evaltmp( nproj), evalcplx( nproj))

  !DO io = 1, nkpt
  !  auxc = zzero
  !  DO ik = 1, nrpt
  !    DO ir = 1, nrpt
  !      auxc = auxc + evaltmp( ik)*exp( zi*DOT_PRODUCT( vkc( :, io)-vkccpy( :, ik), rpt( :, ir)))
  !    END DO
  !  END DO
  !  evalsv( 1, io) = Real( Real( auxc/nrpt))-efermi
  !END DO

  DO ik = 1, nkpt 
    !write(*,*) vkl( :, ik)
    auxmat(:,:) = zzero
    DO ir = 1, nrpt
      auxmat(:,:) = auxmat(:,:) + wanme( ir, :, :)*exp( zi*DOT_PRODUCT( rpt( :, ir), vkc( :, ik)))
    END DO
    Call diaghermat( nproj, auxmat, evaltmp, evectmp)
    evalsv( :, ik) = evaltmp(:)-efermi
    !write(*,*) ik
    !DO io = 1, nproj
    !  DO ia = 1, nproj
    !    write(*,'(F12.5,F12.5," -- ")', advance='no') evectmp( io, ia)
    !  END DO
    !  write(*,*)
    !END DO
    !Call getevalfv( vkl( :, ik), evalfv)
    !DO iproj = 1, nproj
    !  write(*,*) evalfv( pbandstart+iproj-1, 1), evaltmp( iproj)
    !END DO
  END DO
  write(*,*) "end interpolation"
                                                !
    !--------------------------------------------
    
    !============================================
    ! output
        Call xml_OpenFile ("bandstructure.xml", xf, replace=.True., &
             & pretty_print=.True.)
        !
        Call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
             &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
        IF ( .Not. input%properties%bandstructure%character) Then
           Open (50, File='BAND_WANNIER.OUT', Action='WRITE', Form='FORMATTED')

           Call xml_NewElement (xf, "bandstructure")
           Call xml_NewElement (xf, "title")
           Call xml_AddCharacters (xf, trim(input%title))
           Call xml_endElement (xf, "title")
           Do ist = 1, nproj
              Call xml_NewElement (xf, "band")
              Do ik = 1, nkpt
                 io = 0
                 DO ia = 1, nrpt
                   IF( norm2( vklcpy( :, ia)-vkl( :, ik)) .LT. 0.01) THEN
                     io = ia
                     exit
                   END IF
                 END DO
                 IF( io.GT.0) THEN
                    Write (50, '(3G18.10)') dpp1d (ik), evalsv (ist, ik), cf( io, pbandstart+ist-1)-efermi
                 ELSE
                    Write (50, '(2G18.10)') dpp1d (ik), evalsv (ist, ik)
                 END IF
                 Call xml_NewElement (xf, "point")
                 Write (buffer, '(5G18.10)') dpp1d (ik)
                 Call xml_AddAttribute (xf, "distance", &
                      & trim(adjustl(buffer)))
                 Write (buffer, '(5G18.10)') evalsv (ist, ik)
                 Call xml_AddAttribute (xf, "eval", &
                      & trim(adjustl(buffer)))
                 Call xml_endElement (xf, "point")
              End Do
              Call xml_endElement (xf, "band")
              Write (50, '("     ")')
           End Do
           Close (50)
           Write (*,*)
           Write (*, '("Info(bandstr):")')
           Write (*, '(" band structure plot written to BAND_WANNIER.OUT")')
        END IF
                                                !
    !--------------------------------------------
   
  !------------------------------!
  ! End of Wannier interpolation !
  !------------------------------!
stop
#ifdef MPI
!  Call MPI_barrier (MPI_COMM_WORLD, ierr)
!  splittfile = .True.
write (*,*) 'wannier+mpi not implemented'
write(*,*) 'stopping now'
stop
  Do ik = firstk (rank), lastk (rank)
#else
    splittfile = .False.
     Do ik = 1, nkpt
#endif
        Allocate (evalfv(nstfv, nspnfv))
        Allocate (evecfv(nmatmax, nstfv, nspnfv))
        Allocate (evecsv(nstsv, nstsv))
        ! solve the first- and second-variational secular equations

        Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
        Call getevecsv (vkl(:, ik), evecsv)
        Call getevalfv (vkl(:, ik), evalfv)

        Deallocate (evalfv, evecfv, evecsv)
        ! end loop over k-points
     End Do

  ! compute the overlap radial integrals
  !  Call olprad
  ! compute the Hamiltonian radial integrals
  !  Call hmlint
  ! compute "relativistic mass"
  !  Call genmeffig
  !  emin = 1.d5
  !  emax = - 1.d5


else
  !------------------------------------------      
  ! Calculate bandstructure by interpolation
  !------------------------------------------
  hybcheck = .false.
  if (associated(input%groundstate%Hybrid)) then
    hybcheck = .true.
    if (input%properties%bandstructure%character) then
      if (rank == 0) then 
        write(*,'(a)') "Warning(bandstr): "
        write(*,'(a)') "    Atom-resolved bandstructure for hybrids is not yet implemented!"
        input%properties%bandstructure%character = .false.
      end if
    end if
  else if (associated(input%groundstate%HartreeFock)) then
    hybcheck = .true.
  end if

if (hybcheck) then
    !--------------------
    ! begin Interpolation 
    !--------------------
    if (rank==0) then
      fname = 'EVALHF.OUT'
      inquire(File=fname, Exist=exist)
      if (.not.exist) then
        write(*,*)'ERROR(bandstr.f90): File EVALHF.OUT does not exist!'
        stop
      end if
      inquire(IoLength=Recl) nkpt0, nstsv0
      open(70, File=fname, Action='READ', Form='UNFORMATTED', &
      &    Access='DIRECT', Recl=Recl)
      read(70, Rec=1) nkpt0, nstsv0
      close(70)
      nstsv = min(nstsv,nstsv0)
      allocate(vkl0(3,nkpt0))
      allocate(ehf(nstsv0,nkpt0))
      allocate(e0(nkpt0,nstsv))
      allocate(e1(nkpt,nstsv))
      if (allocated(evalsv)) deallocate(evalsv)
      allocate(evalsv(nstsv,nkpt))
      inquire(IoLength=Recl) nkpt0, nstsv0, vkl0(:,1), ehf(:,1)
      open(70, File=fname, Action='READ', Form='UNFORMATTED', &
      &    Access='DIRECT', Recl=Recl)
      do ik = 1, nkpt0
        read(70, Rec=ik) nkpt0, nstsv0, vkl0(:,ik), ehf(:,ik)
      end do ! ik
      close(70)
      ! read fermi energy
      call readfermi
      ! Perform Fourier Interpolation
      do ik = 1, nkpt0
        e0(ik,1:nstsv) = cmplx(ehf(1:nstsv,ik),0.d0,8)
      end do
      e1(:,:) = zzero
      ! Fourier interpolation
      call fourintp(e0,nkpt0,vkl0,e1,nkpt,vkl,nstsv)
      emin =  1.d5
      emax = -1.d5
      do ist = 1, nstsv
        do ik = 1, nkpt
          evalsv(ist,ik) = dble(e1(ik,ist))
          emin = min(emin, evalsv(ist, ik))
          emax = max(emax, evalsv(ist, ik))
        end do
      end do
      deallocate(vkl0,ehf,e0,e1)
      emax = emax + (emax-emin) * 0.5d0
      emin = emin - (emax-emin) * 0.5d0
    end if
    !--------------------
    ! end Interpolation 
    !--------------------

else

  ! maximum angular momentum for band character
  lmax = Min (3, input%groundstate%lmaxapw)
  lmmax = (lmax+1) ** 2
  If (input%properties%bandstructure%character) Then
     Allocate (bc(0:lmax, natmtot, nstsv, nkpt))
  End If
  ! read density and potentials from file
  Call readstate
  ! read Fermi energy from file
  Call readfermi
  ! find the new linearisation energies
  Call linengy
  ! generate the APW radial functions
  Call genapwfr
  ! generate the local-orbital radial functions
  Call genlofr
  ! compute the overlap radial integrals
  Call olprad
  ! compute the Hamiltonian radial integrals
  Call hmlint
  ! compute "relativistic mass"
  Call genmeffig
  emin = 1.d5
  emax = - 1.d5

  !---------------------------------------
  ! begin parallel loop over k-points
  !---------------------------------------
#ifdef MPI
  Call MPI_barrier (MPI_COMM_WORLD, ierr)
  splittfile = .True.
  Do ik = firstk (rank), lastk (rank)
#else
    splittfile = .False.
     Do ik = 1, nkpt
#endif
        Allocate (evalfv(nstfv, nspnfv))
        Allocate (evecfv(nmatmax, nstfv, nspnfv))
        Allocate (evecsv(nstsv, nstsv))
        ! solve the first- and second-variational secular equations
        Call seceqn (ik, evalfv, evecfv, evecsv)
        Do ist = 1, nstsv
           ! subtract the Fermi energy
           evalsv (ist, ik) = evalsv (ist, ik) - efermi
           ! add scissors correction
           If (evalsv(ist, ik) .Gt. 0.d0) evalsv (ist, ik) = evalsv (ist, ik) + &
                & input%properties%bandstructure%scissor
           emin = Min (emin, evalsv(ist, ik))
           emax = Max (emax, evalsv(ist, ik))
        End Do
        ! compute the band characters if required
        If (input%properties%bandstructure%character) Then
           Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
           Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, &
                & nspnfv))
           ! find the matching coefficients
           Do ispn = 1, nspnfv
              Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, &
                   & ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, &
                   & ispn))
           End Do
           ! average band character over spin and m for all atoms
           Do is = 1, nspecies
              Do ia = 1, natoms (is)
                 ias = idxas (ia, is)
                 ! generate the diagonal of the density matrix
                 Call gendmat (.True., .True., 0, lmax, is, ia, ngk(:, &
                      & ik), apwalm, evecfv, evecsv, lmmax, dmat)
                 Do ist = 1, nstsv
                    Do l = 0, lmax
                       sum = 0.d0
                       Do m = - l, l
                          lm = idxlm (l, m)
                          Do ispn = 1, nspinor
                             sum = sum + dble (dmat(lm, lm, ispn, &
                                  & ispn, ist))
                          End Do
                       End Do
                       bc (l, ias, ist, ik) = real (sum)
                    End Do
                 End Do
              End Do
           End Do
           Deallocate (dmat, apwalm)
        End If
        Deallocate (evalfv, evecfv, evecsv)
        ! end loop over k-points
     End Do

#ifdef MPI
     If (input%properties%bandstructure%character) Then
       Call mpi_allgatherv_ifc(nkpt,(lmax+1)*natmtot*nstsv,rlpbuf=bc)
     End If
     Call mpi_allgatherv_ifc(nkpt,nstsv,rbuf=evalsv)
     Call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif

      if (allocated(meffig)) deallocate(meffig)
      if (allocated(m2effig)) deallocate(m2effig)
     emax = emax + (emax-emin) * 0.5d0
     emin = emin - (emax-emin) * 0.5d0
     end if
     if (rank==0) then
        ! output the band structure
        Call xml_OpenFile ("bandstructure.xml", xf, replace=.True., &
             & pretty_print=.True.)
        !
        Call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
             &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
        If ( .Not. input%properties%bandstructure%character) Then
           Open (50, File='BAND.OUT', Action='WRITE', Form='FORMATTED')

           Call xml_NewElement (xf, "bandstructure")
           Call xml_NewElement (xf, "title")
           Call xml_AddCharacters (xf, trim(input%title))
           Call xml_endElement (xf, "title")
           Do ist = 1, nstsv
              Call xml_NewElement (xf, "band")
              Do ik = 1, nkpt
                 Write (50, '(2G18.10)') dpp1d (ik), evalsv (ist, ik)
                 Call xml_NewElement (xf, "point")
                 Write (buffer, '(5G18.10)') dpp1d (ik)
                 Call xml_AddAttribute (xf, "distance", &
                      & trim(adjustl(buffer)))
                 Write (buffer, '(5G18.10)') evalsv (ist, ik)
                 Call xml_AddAttribute (xf, "eval", &
                      & trim(adjustl(buffer)))
                 Call xml_endElement (xf, "point")
              End Do
              Call xml_endElement (xf, "band")
              Write (50, '("     ")')
           End Do
           Close (50)
           Write (*,*)
           Write (*, '("Info(bandstr):")')
           Write (*, '(" band structure plot written to BAND.OUT")')
        Else
           Call xml_NewElement (xf, "bandstructure")
           Call xml_AddAttribute (xf, "character", "true")
           Call xml_NewElement (xf, "title")
           Call xml_AddCharacters (xf, trim(input%title))
           Call xml_endElement (xf, "title")
           Do is = 1, nspecies
              Call xml_NewElement (xf, "species")
              Call xml_AddAttribute (xf, "name", trim(spname(is)))
              Call xml_AddAttribute (xf, "chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol))
              Do ia = 1, natoms (is)
                 Call xml_NewElement (xf, "atom")
                 Write (buffer, '(5G18.10)') atposc (:, ia, is)
                 Call xml_AddAttribute (xf, "coord", &
                      & trim(adjustl(buffer)))
                 ias = idxas (ia, is)
                 Write (fname, '("BAND_S", I2.2, "_A", I4.4, ".OUT")') &
                      & is, ia
                 Open (50, File=trim(fname), Action='WRITE', Form='FORMAT&
                      &TED')
                 !
                 Do ist = 1, nstsv
                    Call xml_NewElement (xf, "band")
                    Do ik = 1, nkpt
                       ! sum band character over l
                       sum = 0.d0
                       Do l = 0, lmax
                          sum = sum + bc (l, ias, ist, ik)
                       End Do
                       Call xml_NewElement (xf, "point")
                       Write (buffer, '(5G18.10)') dpp1d (ik)
                       Call xml_AddAttribute (xf, "distance", &
                            & trim(adjustl(buffer)))
                       Write (buffer, '(5G18.10)') evalsv (ist, ik)
                       Call xml_AddAttribute (xf, "eval", &
                            & trim(adjustl(buffer)))
                       Write (buffer, '(5G18.10)') sum
                       Call xml_AddAttribute (xf, "sum", &
                            & trim(adjustl(buffer)))
                       Do l = 0, lmax
                          Call xml_NewElement (xf, "bc")
                          Write (buffer,*) l
                          Call xml_AddAttribute (xf, "l", &
                               & trim(adjustl(buffer)))
                          Write (buffer, '(5G18.10)') bc (l, ias, ist, &
                               & ik)
                          Call xml_AddAttribute (xf, "character", &
                               & trim(adjustl(buffer)))
                          Call xml_endElement (xf, "bc")
                       End Do
                       Call xml_endElement (xf, "point")
                       Write (50, '(2G18.10, 8F12.6)') dpp1d (ik), evalsv &
                            & (ist, ik), sum, (bc(l, ias, ist, ik), l=0, lmax)
                    End Do
                    Call xml_endElement (xf, "band")
                    Write (50, '("	  ")')
                 End Do
                 Call xml_endElement (xf, "atom")
                 Close (50)
              End Do
              Call xml_endElement (xf, "species")
           End Do
           Write (*,*)
           Write (*, '("Info(bandstr):")')
           Write (*, '(" band structure plot written to BAND_Sss_Aaaaa.OU&
                &T")')
           Write (*, '("	for all species and atoms")')
        End If
        Write (*,*)
        Write (*, '(" Fermi energy is at zero in plot")')
        ! output the vertex location lines
        Open (50, File='BANDLINES.OUT', Action='WRITE', Form='FORMATTED')
        Do iv = 1, nvp1d
           Call xml_NewElement (xf, "vertex")
           !
           Write (buffer, '(5G18.10)') dvp1d (iv)
           Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
           Write (buffer, '(5G18.10)') emax
           Call xml_AddAttribute (xf, "upperboundary", &
                & trim(adjustl(buffer)))
           Write (buffer, '(5G18.10)') emin
           Call xml_AddAttribute (xf, "lowerboundary", &
                & trim(adjustl(buffer)))
           Call xml_AddAttribute (xf, "label", trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label)))
           Write (buffer, '(5G18.10)') input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
           Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
           Call xml_endElement (xf, "vertex")
           Write (50, '(2G18.10)') dvp1d (iv), emin
           Write (50, '(2G18.10)') dvp1d (iv), emax
           Write (50, '("     ")')
        End Do
        Close (50)
        Write (*,*)
        Write (*, '(" vertex location lines written to BANDLINES.OUT")')
        Write (*,*)
        Call xml_close (xf)
     End if
     If (input%properties%bandstructure%character) deallocate(bc)

hybcheck=.false.
if (associated(input%groundstate%Hybrid)) then
    if (input%groundstate%Hybrid%exchangetypenumber== 2)   hybcheck=.true.
elseif (input%groundstate%xctypenumber .Lt. 0) then
       hybcheck=.true.
end if
if (hybcheck) then

     !----------------------------------------      
     ! Calculate interpolated DELTAX
     !----------------------------------------

      
        fname='DELTAX.OUT'
        inquire(File=fname,Exist=exist)
        if (.not.exist) then
          write(*,*)'ERROR(bandstr.f90): File DELTAX.OUT does not exist!'
          stop
        end if

        open(500,file='DELTAX.OUT',action='READ',form='FORMATTED')
        read(500,*) nkpt0, nstsv0

        nstsv=min(nstsv,nstsv0)
        allocate(vkl0(3,nkpt0))
        allocate(deltax(nkpt0,nstsv0))
        allocate(e0(nkpt0,nstsv))
        allocate(e1(nkpt,nstsv))

        do ik=1,nkpt0
              read(500,*) vkl0(1,ik),vkl0(2,ik),vkl0(3,ik)
              read(500,*) deltax(ik,1:nstsv0)
        end do
        close(500)
  
        ! Perform Fourier Interpolation
        do ik = 1, nkpt0
           e0(ik,1:nstsv)=cmplx(deltax(ik,1:nstsv0),0.d0,8)
        enddo

        e1(:,:)=zzero
        call fourintp(e0,nkpt0,vkl0,e1,nkpt,vkl,nstsv)
        
        open(87,file='DELTAX_INTP.OUT')
        do ist = 1, nstsv
          do ik = 1, nkpt
            write(87,*) dpp1d(ik), dble(e1(ik,ist))
          end do
          write(87,*)
        end do
        close(87)

        ! add discontinuity times ex_coef to bandstructure
        open(89,file='BAND_DELTAX.OUT')
        do ist = 1, nstsv
          do ik = 1, nkpt
            write(89,*) dpp1d(ik), evalsv (ist, ik)+ex_coef*dble(e1(ik,ist))
          end do
          write(89,*)
        end do
        close(89)
        deallocate(vkl0,deltax,e0,e1)           
        return
      end if 


endif


    !---------------------------------------------------------------------------
    ! din: New output file for the bandstructure to be able to post-process it
    !---------------------------------------------------------------------------
    if (rank==0) then
      open(50, File="bandstructure.dat", Action='Write', Form='Formatted')
      write(50,*) "# ", 1, nstsv, nkpt
      ! path, energy, ist, ik, vkl
      do ist = 1, nstsv
      do ik = 1, nkpt
        write(50,'(2I6, 3F12.6, 2G18.10)') ist, ik, vkl(:,ik), dpp1d(ik), evalsv(ist,ik)
      end do
      write(50,*)
      end do
      close(50)
    end if

    Return
End Subroutine bandstr
!EOC
