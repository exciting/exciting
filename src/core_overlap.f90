! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bandstr
! !INTERFACE:
!
!
Subroutine core_overlap
  ! !USES:
  Use modinput, only: input
  Use mod_eigensystem, only: mt_hscf, nmatmax, MTNullify, MTInitall 
  Use mod_kpoint, only: nkpt
  Use mod_Gkvector, only: ngkmax, ngk, gkc, tpgkc, sfacgk
  Use mod_APW_LO, only: apwordmax
  Use mod_eigenvalue_occupancy, only: nstfv, nstsv
  Use constants, only: zzero
  Use mod_atoms, only: spr, spl, natmtot, idxas
  Use mod_muffin_tin, only: nrcmtmax, lmmaxapw, nrmt, idxlm, nrcmt, rcmt
  Use modmpi, only: rank, barrier, ierr, firstk, lastk, mpi_allgatherv_ifc 
#ifdef MPI
  Use mpi, only: MPI_COMM_WORLD, MPI_barrier 
#endif
  Use FoX_wxml
  Use mod_spin, only: nspnfv
  Use modxas, only: preml, ucore, ecore, mj2ml, spj, mj, ncg
  Use mod_hdf5

  ! !DESCRIPTION:
  ! Calculates the overlap between a selected group of core states and the 
  ! valence and conduction states. While we assume that the overlap between 
  ! these states vanishes, there is finite overlap since the two groups of 
  ! states are obtained from different Hamiltonians. The overlap 
  ! $\int dr^3 \psi^*_{\mu}(\mathbf{r})\psi_{i\mathbf{k}}(\mathbf{r})$ between 
  ! a core state $\psi_{\mu \mathbf{k}}$ and valence/conduction state 
  ! $\psi_{i\mathbf{k}}$ and energy differences $\Delta \epsilon=\epsilon_{i
  ! \mathbf{k}}-\epsilon_{\mu}$ are stored either in xml or hdf5 file.
  ! !REVISION HISTORY:
  !   Created December 2020 (Christian Vorwerk) based on bandstr.f90
  !EOP
  !BOC
  Implicit None
  ! local variables
  Integer :: lxas
  Integer :: ik, is, ia, ias, ist1, ist2, irc, ir
  Integer :: lm1, lm2
  Real (8), Allocatable :: evalfv (:,:)
  Complex (8), Allocatable :: apwalmt (:, :, :, :)
  Complex (8), Allocatable :: evecfv (:, :)
  Complex (8), Allocatable :: evecsv (:, :)
  Complex (8), Allocatable :: wfmt(:,:,:)
  ! core-valence overlap array
  Complex (8), Allocatable :: overlap(:,:,:)
  ! core-valence energy difference 
  Real(8), Allocatable     :: de(:,:,:)
  Real(8), Allocatable :: fr0(:), fr1(:), fr2(:), gr(:)
  Real(8), Allocatable :: cf(:,:)
  Real(8) :: t1, t2
  character(*), parameter :: thisname="writeoverlap"
  character(256) :: cik
  Character (128) :: buffer
  Type (xmlf_t), Save :: xf

  ! initialise global variables
  Call init0
  ! k-point setup
  Call init1

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
  ! initialize core states
  call coreinit(input%properties%coreoverlap%corespecies, &
    & input%properties%coreoverlap%coreatom)
  
  ! compute "relativistic mass"
  Call genmeffig
  
  allocate(apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
  allocate(wfmt(lmmaxapw,nrcmtmax,nstfv))
  allocate(overlap(ncg,nstfv,nkpt))
  allocate(de(ncg,nstfv,nkpt))
  overlap(:,:,:)=zzero
  de(:,:,:)=0.0d0
  apwalmt(:, :, :, :) = zzero

  ! allocate radial integrals
  Allocate(fr0(nrcmtmax))
  Allocate(fr1(nrcmtmax))
  Allocate(fr2(nrcmtmax))
  Allocate(gr(nrcmtmax))
  Allocate(cf(3,nrcmtmax))
  
  ! allocate eigenstates and -energies  
  Allocate (evalfv(nstfv,1))
  Allocate (evecfv(nmatmax, nstfv))
  Allocate (evecsv(nstsv, nstsv))

  ia=input%properties%coreoverlap%coreatom
  is=input%properties%coreoverlap%corespecies
  ias=idxas(ia,is)
  !---------------------------------------!
  ! begin parallel loop over k-points     !
  !---------------------------------------!
#ifdef MPI
!TODO(Alex) Issue #23 Refactor mpi module
  Do ik = firstk(rank, nkpt), lastk(rank, nkpt)
#else
  Do ik = 1, nkpt
#endif

    ! initialise the eigenvectors if we use the Davidson eigensolver
    if (input%groundstate%solver%type.eq.'Davidson') evecfv=zzero

    ! solve the first- and second-variational secular equations
    Call seceqn (ik, evalfv, evecfv, evecsv)

    ! find the matching coefficients
    Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
     &  sfacgk(1, 1, 1, ik), apwalmt)
    
    ! Calculate the overlap between the core and conduction and valence states
    Do ist1= 1, nstfv
     
     ! calculate the valence or conduction wavefunction
     call wavefmt(input%groundstate%lradstep, &
      &  input%groundstate%lmaxapw,is,ia,ngk(1,ik),apwalmt, &
      &  evecfv(:,ist1),lmmaxapw,wfmt(:,:,ist1))
      Do ist2= 1, ncg
        
        de(ist2, ist1, ik)=evalfv(ist1,1)-ecore(ist2) 
        lxas=spl(ist2,is)
        irc=0
        Do ir=1,nrmt(is),input%groundstate%lradstep
        !Do ir=1,nrmt(is)
          irc=irc+1
          fr0(irc)=ucore(ir,ist2)*spr(ir,is)*spr(ir,is)
        Enddo ! ir
        lm1=idxlm(lxas,mj2ml(lxas,mj(ist2),1))
        lm2=idxlm(lxas,mj2ml(lxas,mj(ist2),2))
        Do irc=1,nrcmt(is)
          fr1(irc)=1.0d0/sqrt(2.0d0)*(preml(lxas,spj(ist2),mj(ist2),1)*fr0(irc)& 
            &   *dble(wfmt(lm1,irc,ist1))+ preml(lxas,spj(ist2),mj(ist2),2)&
            &   *fr0(irc)*dble(wfmt(lm2,irc,ist1)))
           fr2(irc)=1.0d0/sqrt(2.0d0)*(preml(lxas,spj(ist2),mj(ist2),1)*fr0(irc)&
             &   *aimag(wfmt(lm1,irc,ist1))+ preml(lxas,spj(ist2),mj(ist2),2)&
             &   *fr0(irc)*aimag(wfmt(lm2,irc,ist1)))
        Enddo ! irc  
        call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
        t1=gr(nrcmt(is))
        call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
        t2=gr(nrcmt(is))
        overlap(ist2,ist1,ik)=cmplx(t1,t2,8)
      Enddo ! ist2
    Enddo !ist1
    ! end loop over k-points
  End Do

  Deallocate (evalfv, evecfv, evecsv)
  ! gather results when using MPI
#ifdef MPI
  Call mpi_allgatherv_ifc(nkpt, nstfv*ncg, rbuf=de)
  Call mpi_allgatherv_ifc(nkpt, nstfv*ncg, zbuf=overlap)
#endif

!TODO(Bene #159): Replace HDF5 legacy wrappers with the new ones.
!#ifdef _HDF5_
!-----------------------------------------------------------------------------!
!                         Output to HDF5 File                                 !
!                                                                             ! 
!-----------------------------------------------------------------------------!
!  if (rank == 0) then
!    if (.not. hdf5_exist_group(fhdf5, "/", "overlap")) then
!      call hdf5_create_group(fhdf5,"/", "overlap")
!    end if
!    if (.not. hdf5_exist_group(fhdf5, "/", "energy-diff")) then
!      call hdf5_create_group(fhdf5,"/", "energy-diff")
!    end if
!    do ik=1,nkpt 
!      write(cik, '(I8.8)') ik
!      ! Write hdf5
!      call hdf5_write(fhdf5,'overlap',cik, overlap(1,1,ik), shape(overlap(:,:,ik)))
!      call hdf5_write(fhdf5,'energy-diff',cik, de(1,1,ik), shape(de(:,:,ik)))
!    end do
!  end if
!#else
!  call terminate('coreoverlap: HDF5 output can not be created because exciting was not compiled with HDF5.')
!#endif 
!  else 
!-----------------------------------------------------------------------------!
!                          Output to XML File                                 !
!                                                                             ! 
!-----------------------------------------------------------------------------!
  if (rank==0) then

    Call xml_OpenFile ("coreoverlap.xml", xf, replace=.True., pretty_print=.True.)
    Call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
                      &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
    Call xml_NewElement (xf, "coreoverlap")
    Write (buffer, '(I8)') nkpt
    Call xml_AddAttribute (xf, "nkpt", trim(adjustl(buffer)))
    Write (buffer, '(I8)') nstfv
    Call xml_AddAttribute (xf, "nstfv", trim(adjustl(buffer)))
    Write (buffer, '(I8)') ncg
    Call xml_AddAttribute (xf, "ncg", trim(adjustl(buffer)))
    
    do ik = 1, nkpt
        Call xml_NewElement (xf, "kpoint")
        ! write k-point index as attribute
        Write (buffer, '(I8)') ik 
        Call xml_AddAttribute (xf, "index", &
                & trim(adjustl(buffer)))
        ! loop over valence and conduction states
        do ist1 = 1, nstfv
          ! loop over core states
          do ist2 = 1, ncg
            Call xml_NewElement (xf, "pair")
            Write (buffer, '(I8)') ist1
            Call xml_AddAttribute (xf, "ist1", trim(adjustl(buffer)))
            Write (buffer, '(I8)') ist2
            Call xml_AddAttribute (xf, "ist2", trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') de (ist2, ist1, ik)
            Call xml_AddAttribute (xf, "de", &
                & trim(adjustl(buffer)))
            Write (buffer, '(2G18.10)') real(overlap(ist2, ist1, ik)), aimag(overlap(ist2, ist1, ik))
            Call xml_AddAttribute (xf, "overlap", trim(adjustl(buffer)))
            Call xml_endElement (xf, "pair")
          enddo
        enddo
        Call xml_endElement (xf, "kpoint")
    enddo
    Call xml_endElement (xf, "coreoverlap")
    Call xml_close (xf)
   endif
!endif
  call mt_hscf%release()
  Deallocate (apwalmt)
  Deallocate (wfmt,overlap)
  Deallocate (fr0, fr1, fr2, gr, cf)

  Return
End Subroutine core_overlap 
!EOC
