!
subroutine dielmat(a,b)
    
    use modmain
    use modinput
    use modmpi
    implicit none
! input parameters    
    integer, intent(IN) :: a, b
! local variables
    integer :: ik, jk, isym
    integer :: ist, jst, iw
    integer :: recl, iostat
    real(8) :: swidth, eji, t1, t2, t3
    real(8) :: v1 (3), v2 (3), v3 (3), sc(3,3)
    complex(8) :: wplas, eta, zt1
    character(256) :: fname
! allocatable arrays
    real(8), allocatable :: w(:)
    real(8), allocatable :: evalsvt(:), occsvt(:)
    complex (8), allocatable :: pmat(:,:,:)
    complex (8), allocatable :: eps(:)
    integer :: wgrid, lspl
    real(8) :: wmax
    integer :: COMM_LEVEL2
    integer :: kpari, kparf, wpari, wparf

! initialise universal variables
    call init0
    call init1
! read Fermi energy from file
    call readfermi

! allocate local arrays
    wgrid = input%properties%dielmat%wgrid
    allocate(w(wgrid))
    allocate(eps(wgrid))

! generate energy grid
    wmax = input%properties%dielmat%wmax
    t1 = wmax/dble(wgrid)
    Do iw = 1, wgrid
        w(iw) = t1*dble(iw-1)
    End Do
    
! smearing factor
    swidth = input%properties%dielmat%swidth
    eta = cmplx(0.d0,swidth)

! find the record length for momentum matrix element file
    allocate(pmat(3,nstsv,nstsv))
    inquire(iolength=recl) pmat
    deallocate(pmat)
    open(50,File='PMAT.OUT',Action='READ',Form='UNFORMATTED',Access='DIRECT', &
      Recl=recl,IOstat=iostat)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(nonlinopt): error opening PMAT.OUT")')
      write(*,*)
      stop
    end if

    kpari = firstofset(rank, nkptnr)
    kparf = lastofset(rank, nkptnr)
#ifdef MPI
    ! create communicator object
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,kpari,rank/nkptnr,COMM_LEVEL2,ierr)
#endif
       
    eps(:) = zzero
    wplas = zzero

! sum over non-reduced k-points
    do ik = kpari, kparf

! equivalent reduced k-point
        call findkpt(vklnr(:,ik), isym, jk)
        
! read momentum matrix elements from direct-access file
        allocate(pmat(3,nstsv,nstsv))
        read(50, Rec=jk) pmat

! rotate the matrix elements from the reduced to non-reduced k-point
        if (isym > 1) then
            sc(:,:)=symlatc(:,:,lsplsymc(isym))
            do ist = 1, nstsv
                do jst = 1, nstsv
                    v1(:)=dble(pmat(:,ist,jst))
                    call r3mv(sc,v1,v2)
                    v1(:)=aimag(pmat(:,ist,jst))
                    call r3mv(sc,v1,v3)
                    pmat(:,ist,jst)=cmplx(v2(:),v3(:),8)
                end do
            end do
        end if        

! read eigenvalues and occupancies from files
        allocate(evalsvt(nstsv),occsvt(nstsv))
        call getevalsv(vkl(:,jk),evalsvt)
        call getoccsv(vkl(:,jk),occsvt)
       
!--------------------------
! INTERBAND CONTRIBUTION
!--------------------------

! sum over valence states
        do ist = 1, nstsv
        if (evalsvt(ist)<efermi) then

! sum over conduction states
            do jst = 1, nstsv
            if (evalsvt(jst)>efermi) then

                eji = evalsvt(jst)-evalsvt(ist)+input%properties%dielmat%scissor
                if (eji<1.0d-8) then
                    zt1 = zzero
                else
                    zt1 = occmax*pmat(a,ist,jst)*conjg(pmat(b,ist,jst))/(eji*eji)
                    do iw = 1, wgrid
                        eps(iw) = eps(iw)-zt1/(w(iw)-eji+eta)+conjg(zt1)/(w(iw)+eji-eta)
                    end do ! iw
                end if

            end if
            end do ! jst

        end if
        end do ! ist

!--------------------------
! INTRABAND CONTRIBUTION
!--------------------------

        if (input%properties%dielmat%intraband) then
            do ist = 1, nstsv
                zt1 = occmax*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                t1 = evalsvt(ist)-efermi
                wplas = wplas+zt1*delta(t1)
            end do 
        end if
            
        deallocate(pmat)
        deallocate(evalsvt,occsvt)
            
    end do ! ik

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, eps, wgrid, &
    &   MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, wplas, 1, &
    &   MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    zt1 = fourpi/(omega*dble(nkptnr))
    eps(:) = zt1*eps(:)
    if (a .eq. b) eps(:) = eps(:)+zone

    call barrier    

! write the dielectric function to file
    if (rank==0) then
        write(fname, '("EPSILON_", 2I1, ".OUT")') a, b
        open(60, file=trim(fname), action='WRITE', form='FORMATTED')
        do iw = 1, wgrid
            write (60, '(2G18.10)') w(iw), dble(eps(iw))
        end do
        Write (60, '("     ")')
        do iw = 1, wgrid
            write (60, '(2G18.10)') w(iw), aimag(eps(iw))
        end do
        close (60)
        write(*,*)
        write(*, '("Info(dielmat):")')
        write(*, '("  Interband dielectric tensor written to ", a)') trim(adjustl(fname))
    end if

    if (input%properties%dielmat%intraband) then

        wpari = firstofset(rank, wgrid)
        wparf = lastofset(rank, wgrid)
#ifdef MPI
        ! create communicator object
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,wpari,rank/wgrid,COMM_LEVEL2,ierr)
#endif        

        wplas = zt1*wplas
! Drude-like shape is assumed
        do iw = wpari, wparf
            t1 = wplas/(w(iw)*w(iw)+swidth*swidth)
            ! Real part
            t2 = 1.d0-t1
            ! Imaginary part
            t3 = swidth*t1/(w(iw)+eta)
            zt1 = cmplx(t2,t3,8)
            eps(iw) = eps(iw)+zt1
        end do

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, eps, wgrid, &
       &  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif        
        
        call barrier
        if (rank==0) then
            ! write the plasma frequency to file
            write(fname, '("PLASMA_", 2I1, ".OUT")') a, b
            write(*, '("  Plasma frequency written to ", a)') trim(adjustl(fname))
            open(60, File=trim(fname), Action='WRITE', Form='FORMATTED')
            write(60, '(2G18.10, " : plasma frequency")') sqrt(wplas)
            close(60)
            ! write the total dielectric function to file
            write(fname, '("EPSILON_TOTAL_", 2I1, ".OUT")') a, b
            write(*, '("  The total dielectric tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            do iw = 1, wgrid
                write(60, '(2G18.10)') w(iw), dble(eps(iw))
            end do
            write(60, '("     ")')
            do iw = 1, wgrid
                write(60, '(2G18.10)') w(iw), aimag(eps(iw))
            end do
            close(60)
            write(*,*)
        end if
     end if
     
     deallocate(w,eps)

contains    

    ! delta function
    real(8) function delta(e)
        implicit none
        real(8), intent(IN) :: e
        real(8) :: t1
        t1 = pi*(e*e+swidth*swidth)
        delta = swidth/t1
    end function delta
    
    ! the Fermi-Dirac function
    real(8) function f0(e,T)
        implicit none
        real(8), intent(IN) :: e, T
        ! Boltzmann constant [Ha/Kelvin]
        real(8), parameter :: kB=3.1668114d-6
        real(8) :: t1
        f0 = 1.d0/(1.d0+dexp(e/(kB*T)))
    end function f0
    
    ! derivative of the Fermi-Dirac function
    real(8) function df0(e,T)
        implicit none
        real(8), intent(IN) :: e, T
        ! Boltzmann constant [Ha/Kelvin]
        real(8), parameter :: kB=3.1668114d-6
        real(8) :: t1, t2
        t1 = 1.0d0/(kB*T)
        t2 = dexp(e*t1)
        df0 = t1*t2/(1.d0+t2)**2
    end function df0
    
end subroutine
