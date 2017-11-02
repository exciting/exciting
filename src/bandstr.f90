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
  use mod_wannier
  use mod_wfint
  use m_wannier_interpolate_density
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
  Integer :: ik, ispn, is, ia, ias, iv, ist, isym
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
  Character (128) :: buffer
  Type (xmlf_t), Save :: xf
  !
  integer :: nkpt0, nstsv0, Recl
  real(8),    Allocatable :: vkl0(:,:), ehf(:,:), deltax(:,:)
  complex(8), allocatable :: e0(:,:), e1(:,:)
  logical :: exist,hybcheck
!WANNIER
  integer :: iq, nv, ix, iy, iz, ngqwf, ngkmaxint, ngkmaxwf, ngqtmp
  real(8) :: s(3), vk(3), v1(3), vl(3), vc(3), dt
  real(8), allocatable :: eval1(:,:), evalint(:,:), dist(:), vvl(:,:)
  real(8), Allocatable :: bc8 (:, :, :, :)
  complex(8), allocatable :: evectmp(:,:), evecint(:,:,:,:)
  type( k_set) :: int_kset
  type( G_set) :: int_Gset
  type( Gk_set) :: int_Gkset
!END WANNIER
! initialise universal variables
  Call init0

  if( input%properties%bandstructure%wannier) then
    !call wfint_interpolate_bandgap
    !stop
    !do ik = 4, 4, 2
    !  !write(*,*) ik
    !  call generate_k_vectors( int_kset, bvec, (/ik, ik, ik/), (/0.d0, 0.d0, 0.d0/), .true.)
    !  !write(*,*) int_kset%nkpt
    !  !call wannier_interpolate_density( int_kset)
    !  call wfint_init( int_kset)
    !  call wfint_interpolate_eigvec
    !  !call wfint_interpolate_occupancy
    !  call wfint_interpolate_density
    !  write(*,'(I3.3,F23.16)') ik, wfint_efermi
    !end do
    !stop
    !call wannier_delfun
    !call wannier_writefun( 1)
    !--------------------------------------------------!      
    ! Calculate bandstructure by Wannier interpolation !
    !--------------------------------------------------!
    if( .not. associated( input%properties%wannier)) then
      write(*,*) " Error (bandstr): Wannier functions have not been calculated."
      call terminate
    end if
  
    write(*,*) 'Interpolate band-structure...'
    
    ! k-points for interpolation
    input%properties%bandstructure%wannier = .false.
    call init1
    call generate_k_vectors( int_kset, bvec, (/1, 1, nkpt/), (/0.d0, 0.d0, 0.d0/), .false.)
    int_kset%vkl = vkl
    int_kset%vkc = vkc
    !call generate_G_vectors( int_Gset, bvec, intgv, input%groundstate%gmaxvr)
    !call generate_Gk_vectors( int_Gkset, int_kset, wf_Gset, gkmax)
    !ngkmaxint = int_Gkset%ngkmax
    !ngkmaxwf = wf_Gkset%ngkmax
    !write(*,'("ngkmax sys = ",I)') ngkmax
    !write(*,'("ngkmax int = ",I)') ngkmaxint
    !write(*,'("ngkmax wan = ",I)') ngkmaxwf
    allocate( evalint( wf_fst:wf_lst, int_kset%nkpt))
    allocate( evecint( nmatmax, nstsv, nspinor, int_kset%nkpt))
    evecint = zzero
    write(*,'("  interpolation grid set up")')

    ! k-points on grid
    input%properties%bandstructure%wannier = .true.
    call init1
    call readkpts
    !select case (input%properties%wannier%input)
    !  case( "groundstate")
    !  case( "gw")
    !    call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
    !    vkl = wf_kset%vkl
    !    vkc = wf_kset%vkc
    !    nkpt = wf_kset%nkpt
    !    call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
    !  case( "hybrid")
    !    call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)
    !    !vkl = wf_kset%vkl
    !    !vkc = wf_kset%vkc
    !    !nkpt = wf_kset%nkpt
    !    !! recalculate G+k-vectors
    !    !do ik = 1, nkpt
    !    !  do is = 1, nspnfv
    !    !    vl (:) = vkl(:, ik)
    !    !    vc (:) = vkc(:, ik)
    !    !    call gengpvec( vl, vc, ngk( is, ik), igkig( :, is, ik), vgkl( :, :, is, ik), vgkc( :, :, is, ik), gkc( :, is, ik), tpgkc( :, :, is, ik))
    !    !    call gensfacgp( ngk(is, ik), vgkc( :, :, is, ik), ngkmax, sfacgk( :, :, is, ik))
    !    !  end do
    !    !end do
    !    call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
    !  case default
    !    write(*, '(" ERROR (wannier_init): ",a," is not a valid input.")') input%properties%wannier%input
    !    call terminate
    !end select

    allocate( eval1( nstfv, wf_kset%nkpt))
    write(*,'("  do interpolation")')
    call wfint_init( int_kset)
    
    ! read Fermi energy from file
    Call readfermi
      
    ! do interpolation
    lmax = min( 5, input%groundstate%lmaxapw)
    lmmax = (lmax+1)**2
    evalint = wfint_eval
    if( input%properties%bandstructure%character) then
      write(*,*) "character true"
      ! interpolate eigenenergies and eigenvectors
      allocate( bc8( 0:lmax, natmtot, wf_fst:wf_lst, int_kset%nkpt))
      call wfint_interpolate_bandchar( lmax, bc8)
    end if
    evalint = evalint - efermi

    ! output
    ! k-points along path
    input%properties%bandstructure%wannier = .false.
    call init1

    call xml_OpenFile ("bandstructure-wannier.xml", xf, replace=.True., &
         & pretty_print=.True.)
    call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
         &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
      
    allocate( dist( wf_kset%nkpt))
    if( .not. input%properties%bandstructure%character) then
      open( 50, file='BAND_WANNIER.OUT', action='WRITE', form='FORMATTED')
      call xml_NewElement( xf, "bandstructure")
      call xml_NewElement( xf, "title")
      call xml_AddCharacters( xf, trim( input%title))
      call xml_endElement( xf, "title")
      do ist = wf_fst, wf_lst
        call xml_NewElement( xf, "band")
        do iq = 1, int_kset%nkpt
          do ik = 1, wf_kset%nkpt
            dist( ik) = norm2( wf_kset%vkc( :, ik) - int_kset%vkc( :, iq))
          end do
          write( 50, '(3G18.10)') dpp1d (iq), evalint (ist, iq), 1.d0/(1.d0+minval( dist))
          call xml_NewElement( xf, "point")
          write( buffer, '(5G18.10)') dpp1d (iq)
          call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
          write( buffer, '(5G18.10)') evalint( ist, iq)
          call xml_AddAttribute (xf, "eval", trim( adjustl( buffer)))
          call xml_endElement( xf, "point")
        end do
        call xml_endElement( xf, "band")
        write( 50, '("     ")')
      end do
      close(50)
      write(*,*)
      write( *, '("Info (bandstr):")')
      write( *, '(" band structure plot written to BAND_WANNIER.OUT")')
      call xml_endElement( xf, "bandstructure")
      call xml_close( xf)
    else
      call xml_NewElement (xf, "bandstructure")
      call xml_AddAttribute (xf, "character", "true")
      call xml_NewElement (xf, "title")
      call xml_AddCharacters (xf, trim(input%title))
      call xml_endElement (xf, "title")
      do is = 1, nspecies
        call xml_NewElement (xf, "species")
        call xml_AddAttribute (xf, "name", trim(spname(is)))
        call xml_AddAttribute (xf, "chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol))
        do ia = 1, natoms (is)
          call xml_NewElement (xf, "atom")
          write (buffer, '(5G18.10)') atposc (:, ia, is)
          call xml_AddAttribute (xf, "coord", &
               & trim(adjustl(buffer)))
          ias = idxas (ia, is)
          write (fname, '("BAND_WANNIER_S", I2.2, "_A", I4.4, ".OUT")') is, ia
          open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
          !
          do ist = wf_fst, wf_lst
            call xml_NewElement (xf, "band")
            do iq = 1, int_kset%nkpt
              ! sum band character over l
              sum = 0.d0
              do l = 0, lmax
                sum = sum + bc8( l, ias, ist, iq)
              end do
              call xml_NewElement (xf, "point")
              write (buffer, '(5G18.10)') dpp1d (iq)
              call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
              write (buffer, '(5G18.10)') evalint( ist, iq)
              call xml_AddAttribute (xf, "eval", trim(adjustl(buffer)))
              write (buffer, '(5G18.10)') sum
              call xml_AddAttribute (xf, "sum", trim(adjustl(buffer)))
              do l = 0, lmax
                call xml_NewElement (xf, "bc")
                write (buffer,*) l
                call xml_AddAttribute (xf, "l", trim(adjustl(buffer)))
                write (buffer, '(5G18.10)') bc8( l, ias, ist, iq)
                call xml_AddAttribute (xf, "character", trim(adjustl(buffer)))
                call xml_endElement (xf, "bc")
              end do
              call xml_endElement (xf, "point")
              write (50, '(2G18.10, 20F12.6)') dpp1d( iq), evalint( ist, iq), sum, (bc8( l, ias, ist, iq), l=0, lmax)
              !write (50, '(2(G18.10,1x), 8(G12.6,1x))') dpp1d (iq), evalint( ist, iq), sum, (bc( l, ias, ist, iq), l=0, lmax)
            end do
            call xml_endElement (xf, "band")
            write (50, '("	  ")')
          end do
          call xml_endElement (xf, "atom")
          close (50)
        end do
        call xml_endElement (xf, "species")
      end do
      call xml_endElement( xf, "bandstructure")
      call xml_close( xf)
      write (*,*)
      write (*, '("Info(bandstr):")')
      write (*, '(" band structure plot written to BAND_WANNIER_Sss_Aaaaa.OUT")')
      write (*, '("	for all species and atoms")')
    end if  
  
    nv = size( input%properties%bandstructure%plot1d%path%pointarray)
    allocate( vvl( 3, nv))
    
    do ik = 1, nv
      vvl( :, ik) = input%properties%bandstructure%plot1d%path%pointarray( ik)%point%coord
    end do

    ! go back on k-grid
    input%properties%bandstructure%wannier = .true.
    call init1
    call readkpts
    !select case (input%properties%wannier%input)
    !  case( "groundstate")
    !  case( "gw")
    !    call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
    !    vkl = wf_kset%vkl
    !    vkc = wf_kset%vkc
    !    nkpt = wf_kset%nkpt
    !    call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
    !  case( "hybrid")
    !    call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)
    !    vkl = wf_kset%vkl
    !    vkc = wf_kset%vkc
    !    nkpt = wf_kset%nkpt
    !    ! recalculate G+k-vectors
    !    do ik = 1, nkpt
    !      do is = 1, nspnfv
    !        vl (:) = vkl(:, ik)
    !        vc (:) = vkc(:, ik)
    !        call gengpvec( vl, vc, ngk( is, ik), igkig( :, is, ik), vgkl( :, :, is, ik), vgkc( :, :, is, ik), gkc( :, is, ik), tpgkc( :, :, is, ik))
    !        call gensfacgp( ngk(is, ik), vgkc( :, :, is, ik), ngkmax, sfacgk( :, :, is, ik))
    !      end do
    !    end do
    !    call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
    !  case default
    !    write(*, '(" ERROR (wannier_init): ",a," is not a valid input.")') input%properties%wannier%input
    !    call terminate
    !end select
    
    allocate( evalfv( nstfv, nspnfv))
    do ik = 1, wf_kset%nkpt
      call getevalfv( wf_kset%vkl( :, ik), evalfv)
      eval1( :, ik) = evalfv( :, 1)
    end do
    open( 50, file='BANDONGRID.OUT', action='WRITE', form='FORMATTED')
    dt = 0.d0
    do ik = 1, nv-1
      v1 = vvl( :, ik+1) - vvl( :, ik)
      do iz = -wf_kset%ngridk(1), wf_kset%ngridk(1)
        do iy = -wf_kset%ngridk(1), wf_kset%ngridk(1)
          do ix = -wf_kset%ngridk(1), wf_kset%ngridk(1)
            vk(1) = dble( ix)/wf_kset%ngridk(1)
            vk(2) = dble( iy)/wf_kset%ngridk(2)
            vk(3) = dble( iz)/wf_kset%ngridk(3)
            s = -1.d0
            if( abs( v1( 1)) .ge. input%structure%epslat) s(1) = (vk(1) - vvl( 1, ik))/v1( 1)
            if( abs( v1( 2)) .ge. input%structure%epslat) s(2) = (vk(2) - vvl( 2, ik))/v1( 2)
            if( abs( v1( 3)) .ge. input%structure%epslat) s(3) = (vk(3) - vvl( 3, ik))/v1( 3)
            !write( *, '(3F13.6)') s
            if( (s(2) .ge. 0.d0) .and. (s(2) .le. 1.d0)) s(1) = s(2)
            if( (s(3) .ge. 0.d0) .and. (s(3) .le. 1.d0)) s(1) = s(3)
            if( (s(1) .ge. 0.d0) .and. (s(1) .le. 1.d0) .and. &
                (norm2( vvl( :, ik) + s(1)*v1 - vk) .lt. input%structure%epslat)) then
              call r3mv( bvec, vk - vvl( :, ik), s)
              call findkpt( vk, is, iv)
              do is = wf_fst, wf_lst
                write( 50, '(8G18.10,I)') dt + norm2( s), eval1( is, iv)-efermi, vk, vkl( :, iv), iv
              end do
              write( 50, *)
            end if
          end do
        end do
      end do
      call r3mv( bvec, v1, s)
      dt = dt + norm2( s)
    end do
    close( 50)
    !call eqpongrid 
  
    ! go back on BZ-path
    input%properties%bandstructure%wannier = .false.
    call init1
    input%properties%bandstructure%wannier = .true.
    !------------------------------!
    ! End of Wannier interpolation !
    !------------------------------!
  !#ifdef MPI
  !!  Call MPI_barrier (MPI_COMM_WORLD, ierr)
  !!  splittfile = .True.
  !write (*,*) 'wannier+mpi not implemented'
  !write(*,*) 'stopping now'
  !stop
  !  Do ik = firstk (rank), lastk (rank)
  !#else
  !    splittfile = .False.
  !     Do ik = 1, nkpt
  !#endif
  !        Allocate (evalfv(nstfv, nspnfv))
  !        Allocate (evecfv(nmatmax, nstfv, nspnfv))
  !        Allocate (evecsv(nstsv, nstsv))
  !        ! solve the first- and second-variational secular equations
  !
  !        Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
  !        Call getevecsv (vkl(:, ik), evecsv)
  !        Call getevalfv (vkl(:, ik), evalfv)
  !
  !        Deallocate (evalfv, evecfv, evecsv)
  !        ! end loop over k-points
  !     End Do
  !
  !  ! compute the overlap radial integrals
  !  !  Call olprad
  !  ! compute the Hamiltonian radial integrals
  !  !  Call hmlint
  !  ! compute "relativistic mass"
  !  !  Call genmeffig
  !  !  emin = 1.d5
  !  !  emax = - 1.d5
  
  
  else
    call init1
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
           apwalm = zzero
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
    if( input%properties%bandstructure%wannier) then
      write(50,*) "# ", wf_fst, wf_lst, nkpt
      ! path, energy, ist, ik, vkl
      do ist = wf_fst, wf_lst
        do ik = 1, nkpt
          write(50,'(2I6, 3F12.6, 2G18.10)') ist, ik, vkl(:,ik), dpp1d(ik), evalint(ist,ik)
        end do
      write(50,*)
      end do
      deallocate( evalint)
    else
      write(50,*) "# ", 1, nstsv, nkpt
      ! path, energy, ist, ik, vkl
      do ist = 1, nstsv
        do ik = 1, nkpt
          write(50,'(2I6, 3F12.6, 2G18.10)') ist, ik, vkl(:,ik), dpp1d(ik), evalsv(ist,ik)
        end do
      write(50,*)
      end do
    end if
    close(50)
  end if

  Return
End Subroutine bandstr
!EOC
