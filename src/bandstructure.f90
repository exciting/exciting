
!> Module for bandstructure calculations
  
module bandstructure
  use constants, only: zzero, twopi
  use modinput, only: input
  use mod_eigensystem, only: releasesingular
  use precision, only: dp, i32, sp, str_128, str_256

  implicit none

  private

  public :: bandstr, fourintp, bandstr_fourintp

contains
  !>   Produce a band structure along the path in reciprocal-space which connects
  !>   the vertices in the array `vvlp1d`. The band structure is obtained from
  !>   the second-variational eigenvalues and is written to the file `BAND.OUT`
  !>   with the Fermi energy set to zero. If required, band structures are plotted
  !>   to files `BAND\_Sss\_Aaaaa.OUT` for atom `aaaa` of species `ss`,
  !>   which include the band characters for each \(l\) component of that atom in
  !>   columns 4 onwards. Column 3 contains the sum over \(l\) of the characters.
  !>   Vertex location lines are written to `BANDLINES.OUT`.  
  !> REVISION HISTORY
  !> - Created June 2003 (JKD)
  !> - Modified June 2012 (DIN)
  !> - Modified March 2014 (UW)
  !> - Modified June 2018 (SeTi)
  !> - Modified April 2025 (Ronaldo)
subroutine bandstr
  use mod_APW_LO, only: apwordmax
  use mod_atoms, only: atposc, idxas, natoms, natmtot, nspecies, spname
  use mod_Gkvector, only: gkc, ngk, ngkmax, sfacgk, tpgkc
  use mod_eigensystem, only: mt_hscf, MTInitAll, MTNullify, nmatmax, singular
  use mod_eigenvalue_occupancy, only: efermi, evalsv, nstfv, nstsv
  use exciting_mpi, only: xmpi_allgatherv
  use mod_kpoint, only: nkpt, vkl
  use mod_muffin_tin, only: idxlm, lmmaxapw
  use mod_plotting, only: dpp1d, dvp1d, nvp1d
  use mod_potential_and_density, only: meffig, m2effig, xctype
  use mod_spin, only: nspinor, nspnfv
  use modmpi, only: barrier, firstofset, lastofset, mpiglobal, rank, splittfile, terminate
  use FoX_wxml, only: xmlf_t, xml_AddAttribute, xml_AddCharacters, xml_AddXMLPI, xml_Close, &
                      xml_EndElement, xml_NewElement, xml_OpenFile
  use m_write_hdf5, only: hdf5_bandstructure_output
  use precision, only: dp 
  use secular_equation, only: seceqn
  use mod_gen_lo, only: genlofr

  Integer(i32) :: lmax, lmmax, l, m, lm, n_kpts_current_rank
  Integer(i32) :: ik, ispn, is, ia, ias, iv, ist
  Real (dp) :: emin, emax, sum
  Character (len=str_256) :: fname
  Real (dp), Allocatable :: evalfv (:, :)
  ! low precision for band character array saves memory
  Real (sp), Allocatable :: bc (:, :, :, :)
  Complex (dp), Allocatable :: dmat (:, :, :, :, :)
  Complex (dp), Allocatable :: apwalm (:, :, :, :, :)
  Complex (dp), Allocatable :: evecfv (:, :, :)
  Complex (dp), Allocatable :: evecsv (:, :)
  Character (len=str_128) :: buffer
  Type (xmlf_t), Save :: xf

  character(:), allocatable :: label_names
  real(dp), allocatable :: label_coordinates(:, :)

  ! initialise global variables
  Call init0
  Call init1
  if ( trim(input%groundstate%solver%type) == 'Davidson' ) call releasesingular

  !------------------------------------------------------------------
  ! In case of hybrid functionals, one is asked to use Wannier tools
  !------------------------------------------------------------------
  if (xctype(1) >= 400) then
    if (rank==0) then
      write(*,*)
      write(*,*) 'ERROR(bandstr): Please employ the WANNIER tool to interpolate the electronic band structure!'
      write(*,*)
    end if
    call terminate()
  end if

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
  call MTNullify(mt_hscf)
  call MTInitAll(mt_hscf)
  call hmlint(mt_hscf)

  ! compute "relativistic mass"
  Call genmeffig
  emin = 1.d5
  emax = - 1.d5

  !---------------------------------------
  ! begin parallel loop over k-points
  !---------------------------------------
#ifdef MPI
  call barrier
  splittfile = .True.
  Do ik = firstofset(rank, nkpt), lastofset(rank, nkpt)
#else
  splittfile = .False.
  Do ik = 1, nkpt
#endif
    Allocate (evalfv(nstfv, nspnfv))
    Allocate (evecfv(nmatmax, nstfv, nspnfv))
    Allocate (evecsv(nstsv, nstsv))
    ! initialise the eigenvectors if we use the Davidson eigensolver
    ! singular array was built for the GS k grid, and not for 
    ! the bandstr one, so we deallocate them
    if ( trim( input%groundstate%solver%type ) == 'Davidson' ) then
      evecfv = zzero
    end if
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
    End If ! character
    Deallocate (evalfv, evecfv, evecsv)
    ! end loop over k-points
  End Do
  call mt_hscf%release()
  call releasesingular
  n_kpts_current_rank = lastofset(rank, nkpt) - firstofset(rank, nkpt) + 1
  if ( input%properties%bandstructure%character ) &
    call xmpi_allgatherv( mpiglobal, bc, (lmax + 1) * natmtot * nstsv * n_kpts_current_rank )
  call xmpi_allgatherv( mpiglobal, evalsv, nstsv * n_kpts_current_rank )

  if (allocated(meffig)) deallocate(meffig)
  if (allocated(m2effig)) deallocate(m2effig)

  emax = emax + (emax-emin) * 0.5d0
  emin = emin - (emax-emin) * 0.5d0

  allocate(label_coordinates(3, nvp1d))
  label_names = trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(1)%point%label))
  label_coordinates(:, 1) = input%properties%bandstructure%plot1d%path%pointarray(1)%point%coord
  do iv=2, nvp1d
    label_names = label_names // "," // trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label))
    label_coordinates(:, iv) = input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
  end do

  if (input%properties%bandstructure%character) then
    call hdf5_bandstructure_output(mpiglobal, 'properties.h5', '/', evalsv, [emin, emax], dpp1d, label_names, dvp1d, label_coordinates, characters=bc)
  else
    call hdf5_bandstructure_output(mpiglobal, 'properties.h5', '/', evalsv, [emin, emax], dpp1d, label_names, dvp1d, label_coordinates)
  end if  

  !------------------------------
  ! output the band structure
  !------------------------------
  if (rank==0) then
    Call xml_OpenFile ("bandstructure.xml", xf, replace=.True., pretty_print=.True.)
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
          Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
          ias = idxas (ia, is)
          Write (fname, '("BAND_S", I2.2, "_A", I4.4, ".OUT")') is, ia
          Open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
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
                  Write (buffer, '(5G18.10)') bc (l,ias,ist,ik)
                  Call xml_AddAttribute (xf, "character", trim(adjustl(buffer)))
                  Call xml_endElement (xf, "bc")
              End Do
              Call xml_endElement (xf, "point")
              Write (50, '(2G18.10, 8F12.6)') dpp1d (ik), evalsv(ist, ik), sum, (bc(l,ias,ist,ik), l=0,lmax)
            End Do
            Call xml_endElement (xf, "band")
            Write (50, '("	  ")')
          End Do
          Call xml_endElement (xf, "atom")
          Close (50)
        End Do ! ia
        Call xml_endElement (xf, "species")
      End Do ! is
      Write (*,*)
      Write (*, '("Info(bandstr):")')
      Write (*, '(" band structure plot written to BAND_Sss_Aaaaa.OUT")')
      Write (*, '("	for all species and atoms")')
    End If

    Write (*,*)
    Write (*, '(" Fermi energy is at zero in plot")')

    ! output the vertex location lines
    Open (50, File='BANDLINES.OUT', Action='WRITE', Form='FORMATTED')
    Do iv = 1, nvp1d
      Call xml_NewElement (xf, "vertex")
      Write (buffer, '(5G18.10)') dvp1d (iv)
      Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
      Write (buffer, '(5G18.10)') emax
      Call xml_AddAttribute (xf, "upperboundary", trim(adjustl(buffer)))
      Write (buffer, '(5G18.10)') emin
      Call xml_AddAttribute (xf, "lowerboundary", trim(adjustl(buffer)))
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
    Write (*, '(" Vertex location lines written to BANDLINES.OUT")')
    Write (*,*)
    Call xml_close (xf)
  end if ! rank

  If (input%properties%bandstructure%character) deallocate(bc)

  !---------------------------------------------------------------------------------------------------------
  ! Sorry! One more (1000+1, :-) output file for the band structure to be able to apply interpolation on it
  !---------------------------------------------------------------------------------------------------------
  if (rank==0) then
    open(50, File="bandstructure.dat", Action='Write', Form='Formatted')
    write(50,*) "# ", 1, nstsv, nkpt
    do ist = 1, nstsv
      do ik = 1, nkpt
        write(50,'(2I6, 3F12.6, 2G18.10)') ist, ik, vkl(:,ik), dpp1d(ik), evalsv(ist,ik)
      end do
    write(50,*)
    end do
    close(50)
  end if
End Subroutine bandstr

!> This subroutine interpolate function \(f1_n( k )\) defined on the kmesh 1, to kmesh 2
!> using 3D Smooth Fourier transform according to PRB 38, 2721 (1988).
subroutine fourintp(f1, nk1, kvecs1, f2, nk2, kvecs2, nb)
  use constants,    only: zone
  use modmpi,       only: terminate_if_false
  use linear_system_ill_defined_safe, only: ill_defined_safe_solve
  use fouri,        only: nrr, nst, rbas, rindex, rst, setrindex_done
  use mod_symmetry, only: nsymcrys
  use xlapack,      only: matrix_multiply

  integer(i32), intent(in) :: nk1,nk2,nb
  real(dp),    intent(in) :: kvecs1(3,nk1),kvecs2(3,nk2)
  complex(dp), intent(in) :: f1(nk1,1:nb)
  complex(dp), intent(out):: f2(nk2,1:nb) 
    
  integer(i32) :: i, ist, ib, ik, jk, ir
  integer(i32) :: info, rank, lwork, lrwork, liwork
  real(dp), parameter   :: zero_tolerance = 1.0e-6_dp
  real(dp), allocatable :: singular_values(:)
  complex(dp), allocatable :: work(:)
  real(dp), allocatable :: rwork(:)
  integer(i32), allocatable :: iwork(:)
    
  real(dp) :: den, pref, kdotr
  real(dp) :: rmin,rlen,x2,x6,c1,c2
  real(dp) :: r(3), rvec(3), kvec(3)
  real(dp), allocatable  :: rho(:)
        
  complex(dp) :: expkr
  complex(dp), allocatable :: dele(:,:)
  complex(dp), allocatable :: h(:,:)
  complex(dp), allocatable :: coef(:,:)
  complex(dp), allocatable :: smat1(:,:),smat2(:,:)
  complex(dp), allocatable :: sm2(:,:)

  logical :: symmetry
  logical, parameter :: symmetry_default = .true.

  ! shortcut for basis vectors 
  rbas(:,1) = input%structure%crystal%basevect(:,1)
  rbas(:,2) = input%structure%crystal%basevect(:,2)
  rbas(:,3) = input%structure%crystal%basevect(:,3)
  
  if( associated(input%gw) ) then
    symmetry = input%gw%symmetryBandstructure
  else
    symmetry = symmetry_default
  end if
  ! N.B. (Ronaldo) Below is an old comment (maybe inaccurate)
  ! disable symmetry (bug somewhere)
  if( .not. symmetry ) nsymcrys = 1

  ! Set rindex
  if (.not.setrindex_done) then  
    call setrindex
    setrindex_done = .true.
  endif
  
  ! roughness coefficients
  c1 = 0.25_dp
  c2 = 0.25_dp

  allocate(smat1(nk1,nst), &
           smat2(nk2,nst), &
           rho(nst),       &
           coef(nst,nb),   &
           singular_values(1:nk1-1),  &
           sm2(1:nk1-1,1:nst), &
           h(1:nk1-1,1:nk1-1), &
           dele(1:nk1-1,1:nb))
  
  den = real(nsymcrys, kind=dp)

  ! Calculate the star expansion function at each irreducible k-point
  smat1(1:nk1,1:nst) = zzero
  do ik = 1, nk1
    kvec(1:3) = kvecs1(1:3,ik)
    do ir = 2, nrr
      ist = rst(1,ir)
      pref = real(rst(2,ir), kind=dp)
      r(1:3) = real(rindex(1:3,ir), kind=dp)
      kdotr = twopi*sum( r(1:3)*kvec(1:3) ) 
      expkr = cmplx( cos(kdotr), sin(kdotr), kind=dp )
      smat1(ik,ist) = smat1(ik,ist)+pref*expkr/den
    enddo
  enddo 

  ! Carefully check the matrix smat1:
  ! it may occasionally have a reduced rank due to linearly dependent rows.
  ! In such cases, the corresponding rows (and thus some k-points) must be discarded.
  
  ! Calculate the curvature function (rho) for each star
  rho(1:nst) = 0.0_dp
  ist = 1
  do ir = 2, nrr
    if (rst(1,ir) /= ist) then
      ist = rst(1,ir)
      r(1:3) = real(rindex(1:3,ir), kind=dp)
      do i = 1, 3
        rvec(i) = r(1)*rbas(i,1)+r(2)*rbas(i,2)+r(3)*rbas(i,3)
      enddo
      rlen = sum(rvec(1:3)*rvec(1:3))
      if (ist == 2) rmin = rlen
      x2 = rlen/rmin
      x6 = x2*x2*x2
      rho(ist) = (1.0_dp-c1*x2)*(1.0_dp-c1*x2)+c2*x6
    endif
  enddo
  
  ! Set sm2(k)=smat(k)-smat(k_nkp) and dele
  do ik = 1, nk1-1
    do ist = 2, nst
      sm2(ik,ist) = smat1(ik,ist)-smat1(nk1,ist)
    enddo
    do ib = 1, nb
      dele(ik,ib) = f1(ik,ib)-f1(nk1,ib)
    enddo
  enddo
  
  ! Calculate the matrix H      
  h(1:nk1-1,1:nk1-1) = zzero
  do ik = 1, nk1-1
    do jk = 1, nk1-1
      do ist = 2, nst
        h(ik,jk) = h(ik,jk)+sm2(ik,ist)*conjg(sm2(jk,ist))/rho(ist)
      enddo
    enddo
  enddo
  
  ! Solve the Linear equations for the Lagrange multipliers
  ! We use SVD to properly manage ill-conditioned H
  call ill_defined_safe_solve(h, dele, threshold=zero_tolerance)

  ! Calculate the coefficients of the Star expansion
  coef(1,1:nb) = f1(nk1,1:nb)
  do ist = 2, nst
    coef(ist,1:nb) = zzero
    do ik = 1, nk1-1
      coef(ist,1:nb) = coef(ist,1:nb)+dele(ik,1:nb)*conjg(sm2(ik,ist))
    enddo
    coef(ist,1:nb) = coef(ist,1:nb)/rho(ist)
    coef(1,1:nb) = coef(1,1:nb)-coef(ist,1:nb)*smat1(nk1,ist)
  enddo
  
  ! Perform an interpolation to the new k-mesh
  smat2(1:nk2,1:nst) = zzero
  do ik = 1, nk2
    kvec(1:3) = kvecs2(1:3,ik)
    do ir = 1, nrr
      ist = rst(1,ir)
      pref = real( rst(2,ir), kind=dp )
      r(1:3) = real( rindex(1:3,ir), kind=dp )
      kdotr = twopi*sum( r(1:3)*kvec(1:3) )
      expkr = cmplx( cos(kdotr), -sin(kdotr), kind=dp )
      smat2(ik, ist) = smat2(ik, ist) + pref*expkr/den
    enddo 
  enddo 
  call matrix_multiply( smat2, coef, f2 )
    
end subroutine

! Equivalent to bandstr but uses
! the smooth Fourier interpolation to compute the
! bandstructure. Alternative to Wannier
! for hybrids and QSGW runs
subroutine bandstr_fourintp
  use mod_eigenvalue_occupancy, only: efermi, nstfv
  use mod_misc,                 only: task
  use mod_kpoint,               only: nkpt, vkl
  use mod_plotting,             only: dpp1d, dvp1d, nvp1d, vplp1d
  use mod_potential_and_density, only: meffig, m2effig
  use mod_spin,                 only: nspnfv
  use modmpi,                   only: rank, terminate, mpiglobal
  use mod_potential_and_density, only: xctype
  use FoX_wxml,                 only: xmlf_t, xml_AddAttribute, xml_AddCharacters, xml_AddXMLPI, &
                                       xml_Close, xml_EndElement, xml_NewElement, xml_OpenFile
  use m_write_hdf5,             only: hdf5_bandstructure_output
  use modinput,                 only: input
  use precision,                only: dp, i32, str_128, str_256

  implicit none

  integer(i32) :: ik, ist, iv, nkpt_bp, nkpt_ks
  real(dp)     :: emin, emax
  real(dp), allocatable :: evalfv(:,:)

  ! Ground-state k-mesh eigenvalues (shifted by Fermi energy)
  real(dp),    allocatable :: kvecs_gs(:, :)
  complex(dp), allocatable :: gs_evalfv(:, :)   ! (nkpt, nstfv) — complex wrapper for fourintp

  ! Band-path k-points and interpolated eigenvalues
  real(dp),    allocatable :: kvecs_bp(:, :)
  complex(dp), allocatable :: evalfv_bp(:, :)   ! (nkpt_bp, nstfv)

  character(len=str_256) :: fname
  character(len=str_128) :: buffer
  character(:), allocatable :: label_names
  real(dp),    allocatable  :: label_coordinates(:, :)

  type(xmlf_t), save :: xf

  ! ----------------------------------------------------------------
  ! 1.  Globals initialization, and allocation of the KS arrays
  ! ----------------------------------------------------------------
  call init0
  ! Here the k-points are the KS ones
  call init1
  call readfermi
  
  nkpt_ks = nkpt
  allocate(kvecs_gs, source = vkl)
  allocate(gs_evalfv(nkpt_ks, nstfv))
  allocate(evalfv(nstfv, nkpt_ks))
  
  emin =  1.0e5_dp
  emax = -1.0e5_dp

  do ik = 1, nkpt_ks
    call getevalfv(kvecs_gs(:,ik), evalfv(:,ik))
    do ist = 1, nstfv
      evalfv(ist, ik) = evalfv(ist, ik) - efermi
      if (evalfv(ist, ik) > 0.0_dp) &
        evalfv(ist, ik) = evalfv(ist, ik) + input%properties%bandstructure%scissor
      emin = min(emin, evalfv(ist, ik))
      emax = max(emax, evalfv(ist, ik))
      ! fourintp expects layout (nk, nb), i.e. k-index first
      gs_evalfv(ik, ist) = cmplx(evalfv(ist, ik), 0.0_dp, kind=dp)
    end do
  end do
  
  if (allocated(meffig))  deallocate(meffig)
  if (allocated(m2effig)) deallocate(m2effig)


  ! Set the task to band. This sets k points globals
  ! to the bandstructure
  task = 20
  call init1

  ! ----------------------------------------------------------------
  ! 2.  Build the band-path k-point list from the plotting mesh.
  ! ----------------------------------------------------------------
  nkpt_bp = size(dpp1d)   ! number of points along the band path

  allocate(kvecs_bp, source = vplp1d)
  allocate(evalfv_bp(nkpt_bp, nstfv))

  ! ----------------------------------------------------------------
  ! 3.  Smooth Fourier interpolation: gs k-mesh → band-path k-mesh
  ! ----------------------------------------------------------------
  call fourintp(gs_evalfv, nkpt_ks, kvecs_gs, evalfv_bp, nkpt_bp, kvecs_bp, nstfv)

  ! ----------------------------------------------------------------
  ! 4.  Collect interpolated eigenvalues and update emin/emax
  ! ----------------------------------------------------------------
  emax = emax + (maxval(real(evalfv_bp,kind=dp)) - minval(real(evalfv_bp,kind=dp))) * 0.5_dp
  emin = emin - (maxval(real(evalfv_bp,kind=dp)) - minval(real(evalfv_bp,kind=dp))) * 0.5_dp

  ! ----------------------------------------------------------------
  ! 5.  Build label arrays (same logic as bandstr)
  ! ----------------------------------------------------------------
  allocate(label_coordinates(3, nvp1d))
  label_names = trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(1)%point%label))
  label_coordinates(:, 1) = input%properties%bandstructure%plot1d%path%pointarray(1)%point%coord
  do iv = 2, nvp1d
    label_names = label_names // "," // &
      trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label))
    label_coordinates(:, iv) = input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
  end do

  ! ----------------------------------------------------------------
  ! 6.  HDF5 output: store the interpolated eigenvalues.
  !     We reallocate the evalfv array (ist, ik layout) expected by the
  !     output routine, so we transpose evalfv_bp back.
  ! ----------------------------------------------------------------
  deallocate(evalfv)
  allocate(evalfv(nstfv,nkpt_bp))
  do ik = 1, nkpt_bp
    do ist = 1, nstfv
      evalfv(ist, ik) = real(evalfv_bp(ik, ist), kind=dp)
    end do
  end do

  call hdf5_bandstructure_output(mpiglobal, 'properties.h5', '/', &
                                  evalfv(:, 1:nkpt_bp), [emin, emax], &
                                  dpp1d, label_names, dvp1d, label_coordinates)

  ! ----------------------------------------------------------------
  ! 7.  Text / XML output  (rank-0 only, mirrors bandstr)
  ! ----------------------------------------------------------------
  if (rank == 0) then

    call xml_OpenFile("bandstructure.xml", xf, replace=.true., pretty_print=.true.)
    call xml_AddXMLPI(xf, "xml-stylesheet", &
      'href="' // trim(input%xsltpath) // &
      '/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
    call xml_NewElement(xf, "bandstructure")
    call xml_AddAttribute(xf, "interpolated", "true")
    call xml_NewElement(xf, "title")
    call xml_AddCharacters(xf, trim(input%title))
    call xml_EndElement(xf, "title")

    open(50, file='BAND.OUT', action='WRITE', form='FORMATTED')

    do ist = 1, nstfv
      call xml_NewElement(xf, "band")
      do ik = 1, nkpt_bp
        write(50, '(2G18.10)') dpp1d(ik), evalfv(ist, ik)
        call xml_NewElement(xf, "point")
        write(buffer, '(5G18.10)') dpp1d(ik)
        call xml_AddAttribute(xf, "distance", trim(adjustl(buffer)))
        write(buffer, '(5G18.10)') evalfv(ist, ik)
        call xml_AddAttribute(xf, "eval", trim(adjustl(buffer)))
        call xml_EndElement(xf, "point")
      end do
      call xml_EndElement(xf, "band")
      write(50, '("     ")')
    end do

    close(50)
    write(*, *)
    write(*, '("Info(bandstr_fourintp):")')
    write(*, '(" Interpolated band structure written to BAND.OUT")')

    ! Vertex location lines
    open(50, file='BANDLINES.OUT', action='WRITE', form='FORMATTED')
    do iv = 1, nvp1d
      call xml_NewElement(xf, "vertex")
      write(buffer, '(5G18.10)') dvp1d(iv)
      call xml_AddAttribute(xf, "distance", trim(adjustl(buffer)))
      write(buffer, '(5G18.10)') emax
      call xml_AddAttribute(xf, "upperboundary", trim(adjustl(buffer)))
      write(buffer, '(5G18.10)') emin
      call xml_AddAttribute(xf, "lowerboundary", trim(adjustl(buffer)))
      call xml_AddAttribute(xf, "label", &
        trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label)))
      write(buffer, '(5G18.10)') &
        input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
      call xml_AddAttribute(xf, "coord", trim(adjustl(buffer)))
      call xml_EndElement(xf, "vertex")
      write(50, '(2G18.10)') dvp1d(iv), emin
      write(50, '(2G18.10)') dvp1d(iv), emax
      write(50, '("     ")')
    end do
    close(50)

    write(*, '(" Vertex location lines written to BANDLINES.OUT")')
    write(*, '(" Fermi energy is at zero in plot")')
    write(*, *)
    call xml_Close(xf)

    ! bandstructure.dat (same format as bandstr, using band-path k indices)
    open(50, file='bandstructure.dat', action='WRITE', form='FORMATTED')
    write(50, *) "# ", 1, nstfv, nkpt_bp
    do ist = 1, nstfv
      do ik = 1, nkpt_bp
        write(50, '(2I6, 3F12.6, 2G18.10)') ist, ik, kvecs_bp(:, ik), dpp1d(ik), evalfv(ist, ik)
      end do
      write(50, *)
    end do
    close(50)
    write(*, '(" Interpolated band structure also written to bandstructure.dat")')

  end if ! rank == 0

  deallocate(evalfv, gs_evalfv, kvecs_gs, kvecs_bp, evalfv_bp, label_coordinates)

end subroutine bandstr_fourintp

end module
