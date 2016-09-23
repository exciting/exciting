
subroutine fermisurf

  use modmain
  use modmpi
  implicit none

! local variables
  integer :: ik,ist,ist0,ist1
  integer :: fnum,fnum0,fnum1
  integer :: lst,nst,i,i1,i2,i3,j1,j2,j3
  real(8) :: vc(3,4)

! allocatable arrays
  real(8), allocatable :: evalfv(:,:)
  complex(8), allocatable :: evecfv(:,:,:)
  complex(8), allocatable :: evecsv(:,:)

! initialise universal variables
  call init0
  call init1
! read density and potentials from file
  call readstate
! read Fermi energy from file
  call readfermi
! find the new linearisation energies
  call linengy
! generate the APW radial functions
  call genapwfr
! generate the local-orbital radial functions
  call genlofr
! compute the overlap radial integrals
  call olprad
! compute the Hamiltonian radial integrals
  call hmlint
! compute "relativistic mass" on the G-grid
  Call genmeffig
! begin parallel loop over k-points
#ifdef MPI
  do ik = firstofset(mod(rank,nkpt),nkpt), lastofset(mod(rank,nkpt),nkpt)
#else
  do ik = 1,nkpt
#endif
    allocate(evalfv(nstfv,nspnfv))
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(evecsv(nstsv,nstsv))
    ! solve the first- and second-variational secular equations
    call seceqn(ik,evalfv,evecfv,evecsv)
    deallocate(evalfv,evecfv,evecsv)
  end do ! ik
#ifdef MPI  
  call mpi_allgatherv_ifc(nkpt,nstsv,rbuf=evalsv)
#endif
  if (allocated(meffig)) deallocate(meffig)
  if (allocated(m2effig)) deallocate(m2effig)

!---------------------------------------------------
! OUTPUT Block  
!---------------------------------------------------
if (rank==0) then

  if (task.Eq.101) then
!---------------------------------------------------
! 2D plot
!---------------------------------------------------
    fnum0=50
    fnum1=50
    if (ndmag.eq.1) then
    ! special case of collinear magnetism
      open(50,file='FERMISURF2D_UP.xsf',action='WRITE',form='FORMATTED')
      open(51,file='FERMISURF2D_DN.xsf',action='WRITE',form='FORMATTED')
      fnum1=51
      ist=nstfv-input%groundstate%nempty
      ist0=max(ist-input%properties%fermisurfaceplot%nstfsp/2,1)
      ist1=min(ist+input%properties%fermisurfaceplot%nstfsp/2,nstfv)
    else
    ! spin-unpolarised and non-collinear cases
      open(50,file='FERMISURF2D.xsf',action='WRITE',form='FORMATTED')
      ist=(nstfv-input%groundstate%nempty)*nspinor
      ist0=max(ist-input%properties%fermisurfaceplot%nstfsp/2,1)
      ist1=min(ist+input%properties%fermisurfaceplot%nstfsp/2,nstsv)
    end if
    nst=ist1-ist0+1
    ! produce the 2D Fermi surface plot           
    do i=1,3
      vc(:,i)=bvec(:,1)*vclp2d(1,i)+bvec(:,2)*vclp2d(2,i)+bvec(:,3)*vclp2d(3,i)
    end do

    lst=0
    do fnum=fnum0,fnum1 
      if ((ndmag.eq.1).and.(fnum.eq.fnum1)) lst=nstfv
      write(fnum,*)'CRYSTAL'
      write(fnum,*)
      write(fnum,*)'PRIMVEC'
      do i = 1, 3
        write(fnum,'(3f14.9)') bvec(:,i) 
      end do
      write(fnum,*)
      write(fnum,*)'PRIMCOORD'
      write(fnum,*)'1 1'
      write(fnum,*)'0 0 0 0'
      write(fnum,*)
      if ((ndmag.eq.1).and.(fnum.eq.fnum1)) lst=nstfv
      do ist=ist0,ist1
        write(fnum,*)'BEGIN_BLOCK_DATAGRID_2D'
        write(fnum,*)'  band_energies #', ist
        write(fnum,*)'BEGIN_DATAGRID_2D'
        write(fnum,'(2i6)') np2d(:)
        do i=1,3
          write(fnum,'(3G18.10)') vc(:,i)
        end do
        do ik = 1, nkpt
          write(fnum,'(7G18.10)') evalsv(ist+lst,ik)-efermi
        end do
        write(fnum,*)'END_DATAGRID_2D'
        write(fnum,*)'END_BLOCK_DATAGRID_2D'
        end do ! ist
      close(fnum)
    end do ! fnum
    write(*,*)
    write(*,'("Info(fermisurf):")')
    if (ndmag.eq.1) then
      write(*,'(" 2D Fermi surface data written to FERMISURF2D_UP.xsf and &
       &FERMISURF2D_DN.xsf")')
    else
      write(*,'(" 2D Fermi surface data written to FERMISURF2D.xsf")')
    end if
    write(*,'(" for plotting with XCrysDen (Fermi energy set to zero)")')
    write(*,*)
    write(*,'(" Launch as: xcrysden --xsf FERMISURF(_UP/_DN).xsf")')

  else
  !--------------------------------------------------
  ! 3D plot
  !--------------------------------------------------
    fnum0=50
    fnum1=50
    if (ndmag.eq.1) then
    ! special case of collinear magnetism
      open(50,file='FERMISURF_UP.bxsf',action='WRITE',form='FORMATTED')
      open(51,file='FERMISURF_DN.bxsf',action='WRITE',form='FORMATTED')
      fnum1=51
      ist=nstfv-input%groundstate%nempty
      ist0=max(ist-input%properties%fermisurfaceplot%nstfsp/2,1)
      ist1=min(ist+input%properties%fermisurfaceplot%nstfsp/2,nstfv)
    else
    ! spin-unpolarised and non-collinear cases
      open(50,file='FERMISURF.bxsf',action='WRITE',form='FORMATTED')
      ist=(nstfv-input%groundstate%nempty)*nspinor
      ist0=max(ist-input%properties%fermisurfaceplot%nstfsp/2,1)
      ist1=min(ist+input%properties%fermisurfaceplot%nstfsp/2,nstsv)
    end if
    nst=ist1-ist0+1
    ! plotting box in Cartesian coordinates
    do i=1,4
      vc(:,i)=bvec(:,1)*vclp3d(1,i)+bvec(:,2)*vclp3d(2,i)+bvec(:,3)*vclp3d(3,i)
    end do
    ! produce the Fermi surface plot
    lst=0
    do fnum=fnum0,fnum1
      write(fnum,'(" BEGIN_INFO")')
      write(fnum,'(" # Band-XCRYSDEN-Structure-File for Fermi surface plotting")')
      write(fnum,'(" # Launch as: xcrysden --bxsf FERMISURF(_UP/_DN).bxsf")')
      write(fnum,'("   Fermi Energy: ",G18.10)') 0.d0
      write(fnum,'(" END_INFO")')
      write(fnum,'(" BEGIN_BLOCK_BANDGRID_3D")')
      write(fnum, '(" band_energies")')
      write(fnum,'(" BANDGRID_3D_BANDS")')
      write(fnum,'(I4)') nst
      write(fnum,'(3I6)') np3d(:)+1
      do i=1,4
        write(fnum,'(3G18.10)') vc(:,i)
      end do
      if ((ndmag.eq.1).and.(fnum.eq.fnum1)) lst=nstfv
      do ist=ist0,ist1
        write(fnum,'(" BAND: ",I4)') ist
        do i1=0,np3d(1)
          do i2=0,np3d(2)
            do i3=0,np3d(3)
              j1=i1; if (i1==np3d(1)) j1=0
              j2=i2; if (i2==np3d(2)) j2=0
              j3=i3; if (i3==np3d(3)) j3=0
              ik=ikmap(j1,j2,j3)
              write(fnum,'(7G18.10)') evalsv(ist+lst,ik)-efermi
            end do
          end do
        end do
      end do
      write(fnum,'(" END_BANDGRID_3D")')
      write(fnum,'(" END_BLOCK_BANDGRID_3D")')
      close(fnum)
    end do
    write(*,*)
    write(*,'("Info(fermisurf):")')
    if (ndmag.eq.1) then
      write(*,'(" 3D Fermi surface data written to FERMISURF_UP.bxsf and &
       &FERMISURF_DN.bxsf")')
    else
      write(*,'(" 3D Fermi surface data written to FERMISURF.bxsf")')
    end if
    write(*,'(" for plotting with XCrysDen (Fermi energy set to zero)")')
    write(*,*)
    write(*,'(" Launch as: xcrysden --bxsf FERMISURF(_UP/_DN).bxsf")')

  end if 

end if ! rank

  return
end subroutine

