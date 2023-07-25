! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl and
! F. Cricchio. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: spintext
! !INTERFACE:
!
!
module spintexture
  ! !USES:
  Use modmain
  Use modmpi
  Use FoX_wxml, only: xmlf_t, xml_OpenFile, xml_NewElement, xml_AddAttribute, xml_EndElement, xml_Close
  use constants, only : zzero

  character(256), parameter :: fname = "spintext.xml"
  contains
  subroutine calculate_spintexture(band_min, band_max)
    ! !DESCRIPTION:
    !   Calculates the spintexture in a grid in the reciprocal space {\tt vclp2d}.
    !   The spintexture is obtained from the second-variational eigenvalues with SOC
    !   and the result is written to the file {\tt spintext.xml}. The spintexture is
    !   defined as the expectation value of the pauli vector at k-points on a grid.
    !   With second variational eigenstates, the expectation value at a given k-point
    !   $\vec{k}$ and state $i$ reads: 
    !   $\vec{s}_i(\vec{k})=\langle \Psi_{i\vec{k}} |\vec{\sigma}| \Psi_{i\vec{k}}\rangle
    !   =\sum_{j=1}^{nstfv} \begin{pmatrix} 2\text{Re}(A^*_{ij}(\vec{k})A_{ij+nstfv}(\vec{k}))\\
    !   -2\text{Im}(A^*_{ij}(\vec{k})A_{ij+nstfv}(\vec{k}))\\|A_{ij}(\vec{k})|^2
    !   -|A_{ij+nftsv}(\vec{k})|^2\end{pmatrix}$
    !   where the second variational wavfunction is given by:
    !   $| \Psi_{i\vec{k}}\rangle=\sum_{l=1}^{nstfv}\{A_{il}(\vec{k})\begin{pmatrix}1\\0
    !   \end{pmatrix} +A_{il+nfstfv}\begin{pmatrix}0\\1\end{pmatrix}\}|\psi_{l\vec{k}}\rangle$
    !   where $|\psi_{l\vec{k}}\rangle$ are the first-variational eigenstates and 
    !   $nstfv$ is the number of these states.
    !   
    ! !REVISION HISTORY:
    !    Created October 2020 (BMM)
    !  
    !EOP
    !BOC
    Implicit None
  ! input/output
    !> minimum and maximum band indices. Calculates spintextures in range band_min to band_max
    Integer, intent(in) :: band_min, band_max 
  ! local variables
    Integer :: ik, j, ist, nst 
    Complex (8)               :: sc1, sc2, sc3, zdotc
    Real    (8), Allocatable  :: stext (: ,: )
    character(256) :: nstsv_str
  ! allocatable arrays
    
    Real    (8), Allocatable  :: evalfv (:, :)
    Complex (8), Allocatable  :: evecsv (:, :), evecfv (:, :, :)
    !
    

  ! initialise universal variables
    Call init0
    Call init1

  ! check if spin is enabled
    if (.not. associated(input%groundstate%spin)) then
      if (rank==0) then
        write(*,*)
        write(*,*)'ERROR(spintexture): Please enable spinorbit coupling in the ground state!'
        write(*,*)
      end if
      call terminate
    end if
    ! check band_min and band_max
    if (band_min<1) then
      if (rank==0) then
        write(*,*)
        write(*,*)'ERROR(spintexture): The lower bound of bands must be at least 1.'
        write(*,*)
      end if
      call terminate
    end if
    if (band_max>nstsv) then
      if (rank==0) then
        write(*,*)
        write (nstsv_str,*) nstsv
        write(*,*)'ERROR(spintexture): The upper bound of bands cannot be greater than nstsv ('//trim(nstsv_str)//').'
        write(*,*)
      end if
      call terminate
    end if

    nst = (band_max-band_min+1)
    Allocate (stext (3*nst, nkpt))
    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,nkpt))
    
    ! read density and potentials from file
    Call readstate 
    ! input: ?
    ! output: rhomt, vclmt, vxcmt, veffmt, rhoir, vclir, vxcir, veffir, veffig
    
    ! read Fermi energy from file
    Call readfermi 
    ! input: ? 
    ! output: efermi
    
    ! find the new linearisation energies
    Call linengy
    ! input: apwe0, lorbe0, veffmt, y00, apwve, nrmt?, spr, efermi
    ! output: apwe, lorbe

    ! generate the APW radial functions
    Call genapwfr
    ! input: veffmt, y00, apwdm, apwe, spr
    ! output: apwfr

    ! generate the local-orbital radial functions
    Call genlofr
    ! input: veffmt, y00, spr
    ! output: lofr

    ! compute the overlap radial integrals
    Call olprad
    ! input: apwfr, lofr
    ! output: h1aa, h1loa, h1lolo

    ! compute the Hamiltonian radial integrals
    call MTNullify(mt_hscf)
    ! input: mt_hscf
    ! output: mt_hscf
    call MTInitAll(mt_hscf)
    ! input: mt_hscf
    ! output: mt_hscf
    call hmlint(mt_hscf)
    ! input: mt_hscf, y00, veffmt, levelzora
    ! output: mt_hscf

    ! compute "relativistic mass"
    call genmeffig
    ! input: veffir, cfunir, ngvec, igfft
    ! output: veffig

      !---------------------------------------
      ! begin parallel loop over k-points
      !---------------------------------------
#ifdef MPI
    Call MPI_barrier (MPI_COMM_WORLD, ierr)
    splittfile = .True.
    Do ik = firstk(rank, nkpt), lastk(rank, nkpt)
#else
    splittfile = .False.
    Do ik = 1, nkpt
#endif
      Allocate (evalfv(nstfv, nspnfv))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (evecsv(nstsv, nstsv))
      ! initialise the eigenvectors if we use the Davidson eigensolver
      if (input%groundstate%solver%type.eq.'Davidson') evecfv=zzero

      ! solve the first- and second-variational secular equations
      Call seceqn (ik, evalfv, evecfv, evecsv)
      ! calculate spin expectation values
      Do ist = band_min, band_max
        sc1 = zdotc(nstfv, evecsv(1:nstfv,ist), 1, evecsv(nstfv+1:nstsv,ist), 1)
        sc2 = zdotc(nstfv, evecsv(1:nstfv,ist), 1, evecsv(1:nstfv,ist), 1)
        sc3 = zdotc(nstfv, evecsv(nstfv+1:nstsv,ist), 1, evecsv(nstfv+1:nstsv,ist), 1)
        
        stext((ist-band_min)*3+1, ik) =  2*dble(sc1)
        stext((ist-band_min)*3+2, ik) = -2*imag(sc1)
        stext((ist-band_min)*3+3, ik) =  dble(sc2-sc3)
      End Do !ist
      deallocate(evalfv, evecfv, evecsv)
    End Do !iki
    call mt_hscf%release()
    
#ifdef MPI
    Call mpi_allgatherv_ifc(nkpt, 3*nst, rbuf=stext)
    Call mpi_allgatherv_ifc(nkpt, nstsv, rbuf=evalsv)
#endif
  If (rank==0) then
  !xml  
    call write_spintexture_xml(fname, stext, evalsv, band_min, band_max)
  End If
      

  End Subroutine calculate_spintexture

  Subroutine write_spintexture_xml(fname, stext, eval, band_min, band_max)
  ! ! DESCRIPTION:
    ! Writes spin texture to an xml file. 
    implicit none
    ! input/output:
    character(256), intent(in) :: fname
    Real(8), intent(in) :: stext (:,:)
    Real(8), intent(in) :: eval(:,:) 
    integer, intent(in) :: band_min, band_max
    ! local variables:
    Integer :: ist, ik, nst, nkpt, shape_stext(2)
    Real :: spinl(3)
    Type (xmlf_t), Save :: xf

    shape_stext = shape(stext)
    nst = shape_stext(1)/3
    nkpt = shape_stext(2)

    Call xml_OpenFile (fname, xf, replace=.True., pretty_print=.True.)
        Call xml_NewElement(xf, "spintext")
        Do ist = 1, nst
          Call xml_NewElement(xf, "band")
            Call xml_AddAttribute(xf, "ist", ist+band_min-1)
          Do ik = 1, nkpt
              spinl(1) = stext((ist-1)*3+1, ik)
              spinl(2) = stext((ist-1)*3+2, ik)
              spinl(3) = stext((ist-1)*3+3, ik)
            Call xml_NewElement(xf, "k-point")
              Call xml_AddAttribute(xf, "energy", (eval(ist+band_min-1, ik)-efermi)*27.2113966413079)
              Call xml_AddAttribute(xf, "vec", vkc (:, ik))
              Call xml_AddAttribute(xf, "spin", spinl (:))
            Call xml_EndElement(xf, "k-point")
          End Do
          Call xml_EndElement(xf, "band")
        End Do 
        Call xml_EndElement(xf, "spintext")
      Call xml_Close(xf)
    End Subroutine write_spintexture_xml

End Module spintexture
