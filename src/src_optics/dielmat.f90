!
subroutine dielmat
    use modmain
    use modinput
    use modmpi
    use modxs, only: symt2
    implicit none
! local variables
    integer :: l, a, b, ncomp, epscomp(2,9)
    integer :: ik, jk, isym
    integer :: ist, jst, iw
    integer :: recl, iostat
    logical :: exist, intraband, drude
    real(8) :: wplas, swidth, eji, t1, t2
    real(8) :: v1 (3), v2 (3), sc(3,3)
    complex(8) :: zt1, zt2
    character(256) :: fname
! allocatable arrays
    real(8), allocatable :: w(:)
    real(8), allocatable :: evalsvt(:), occsvt(:)
    complex (8), allocatable :: pmat(:,:,:)
    complex (8), allocatable :: sigma(:)
    complex (8), allocatable :: eta(:)
    integer :: wgrid
    integer :: kfirst, klast
    real(8) :: wmax, scissor
! external functions
    real(8) sdelta
    
    if (rank==0) then
      write(*,*)
      write(*,'("Info(dielmat):")')
    end if
    
! initialise universal variables
    call init0
    call init1

! read Fermi energy from file
    call readfermi

! working arrays
    wgrid = input%properties%dielmat%wgrid
    allocate(w(wgrid))
    allocate(sigma(wgrid))

! generate energy grid
    wmax = input%properties%dielmat%wmax
    t1 = wmax/dble(wgrid)
    Do iw = 1, wgrid
        w(iw) = t1*dble(iw-1)
    End Do

! if element "drude" present, use the specified parameters for the Drude term 
! instead of calculating the plasma frequency
    drude = .false.
    if (input%properties%dielmat%intraband) then
        ! if omega_p is different from default (negative) value
        if (input%properties%dielmat%drude(1)>1.d-8) then
            drude = .true.
            if (rank==0) then
                write(*,*) 
                write(*,*) '  Drude model parameters:'
                write(*,*) '     plasma frequency   : ', input%properties%dielmat%drude(1)
                write(*,*) '     lifetime broadening: ', input%properties%dielmat%drude(2)
            end if
        end if
    end if

! input%properties%dielmat%scissor
    scissor = input%properties%dielmat%scissor

! smearing factor
    swidth = input%properties%dielmat%swidth

! finite lifetime broadening to mimic the experimental resolution
! chosen according to PRB86, 125139 (2012)
    allocate(eta(wgrid))
    do iw = 1, wgrid
        eta(iw) = cmplx(0.d0,swidth)
        !eta(iw) = cmplx(0.d0,swidth+0.1*swidth*w(iw))
    end do

! calculate momentum matrix elements
    !if (rank==0) then
    !    write(*,*)
    !    write(*,'("  Calculate the momentum matrix elements")')
    !end if
    !call writepmat

! the momentum matrix elements
    allocate(pmat(3,nstsv,nstsv))
    inquire(iolength=recl) pmat
    open(50,File='PMAT.OUT',Action='READ',Form='UNFORMATTED', &
    &    Access='DIRECT',Recl=recl,IOstat=iostat)

! KS eigenvalues and occupancies    
    allocate(evalsvt(nstsv),occsvt(nstsv))

#ifdef MPI
    kfirst = firstofset(rank,nkptnr)
    klast = lastofset(rank,nkptnr)
#else
    kfirst = 1
    klast = nkptnr
#endif

!-------------------------------------------------------------------------------
!   Loop over optical components
!-------------------------------------------------------------------------------

    ncomp = size(input%properties%dielmat%epscomp,2)
    if (ncomp>0) then
        do l = 1, ncomp
            epscomp(1,l) = input%properties%dielmat%epscomp(1,l)
            epscomp(2,l) = input%properties%dielmat%epscomp(2,l)
        end do
    else
!
! calculate \eps_{ij}} for all symmetry inequivalent components,
! by analysing the symmetrization matrices
! 
        ncomp = 0
        do a = 1, 3
            do b = 1, 3
                if (sum(abs(symt2(a,b,:,:)))>1.0d-8) then
                    ncomp = ncomp+1
                    epscomp(1,ncomp) = a
                    epscomp(2,ncomp) = b
                end if
            end do
        end do
    end if

    do l = 1, ncomp
        
        a = epscomp(1,l)
        b = epscomp(2,l)

! flag for calculating the intraband term
        intraband = input%properties%dielmat%intraband .and. (a==b)

! sum over non-reduced k-points
        sigma(:) = zzero
        wplas = 0.0d0
        do ik = kfirst, klast

! equivalent reduced k-point
            call findkpt(vklnr(:,ik),isym,jk)
        
! read momentum matrix elements from direct-access file
            read(50,Rec=jk) pmat

! rotate the matrix elements from the reduced to non-reduced k-point
            if (isym > 1) then
                do ist = 1, nstsv
                    do jst = 1, nstsv
                        sc(:,:)=symlatc(:,:,lsplsymc(isym))
                        call r3mv(sc,dble(pmat(:,ist,jst)),v1)
                        call r3mv(sc,aimag(pmat(:,ist,jst)),v2)
                        pmat(:,ist,jst) = cmplx(v1(:),v2(:),8)
                    end do
                end do
            end if
        
! read eigenvalues and occupancies from files
            call getevalsv(vkl(:,jk),evalsvt)
            call getoccsv(vkl(:,jk),occsvt)

!--------------------------
! PLASMA FREQUENCY
!--------------------------
            if (intraband.and.(.not.drude)) then
                do ist = 1, nstsv
                    zt1 = occsvt(ist)*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                    t1 = (evalsvt(ist)-efermi)/input%groundstate%swidth
                    wplas = wplas+dble(zt1)*sdelta(input%groundstate%stypenumber,t1)/input%groundstate%swidth
                end do
            end if

!--------------------------
! INTERBAND CONTRIBUTION
!--------------------------       
            do ist = 1, nstsv
            if (evalsvt(ist)<=efermi) then
            
                do jst = 1, nstsv
                if (evalsvt(jst)>efermi) then
            
                    zt1 = pmat(a,ist,jst)*conjg(pmat(b,ist,jst))
                    eji = evalsvt(jst)-evalsvt(ist)
                
                    ! scissor operator
                    if (dabs(scissor)>1.0d-8) then
                        if (dabs(eji)>1.0d-8) then
                            t1 = (eji+scissor)/eji
                            zt1 = zt1*t1*t1
                            eji = eji+scissor
                        end if
                    end if
                
                    if (dabs(eji)>1.0d-8) then
                        t1=occsvt(ist)*(1.0d0-occsvt(jst)/occmax)/eji
                        do iw = 1, wgrid
                            sigma(iw) = sigma(iw)+ &
                           &  t1*(zt1/(w(iw)-eji+eta(iw)) + &
                           &  conjg(zt1)/(w(iw)+eji+eta(iw)))
                        end do ! iw
                    end if

                end if
                end do ! jst
        
            end if
            end do ! ist
            
        end do ! ik
        
        zt1 = zi/(omega*dble(nkptnr))
        sigma(:) = zt1*sigma(:)

#ifdef MPI
        call MPI_AllReduce(MPI_IN_PLACE, &
        &                  sigma, &
        &                  wgrid, &
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        if (intraband.and.(.not.drude)) then
          call MPI_AllReduce(MPI_IN_PLACE, &
          &                  wplas, &
          &                  1, &
          &                  MPI_DOUBLE_PRECISION, &
          &                  MPI_SUM, &
          &                  MPI_COMM_WORLD, &
          &                  ierr)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

! add the intraband contribution
        if (intraband) then
            if (drude) then
! use specified parameters for Drude term
                zt1 = zi/fourpi
                do iw = 1, wgrid
                    t1 = input%properties%dielmat%drude(1)
                    t2 = input%properties%dielmat%drude(2)
                    sigma(iw) = sigma(iw)+zt1*t1*t1/(w(iw)+zi*t2)
                end do
            else
! use the calculated plasma frequency            
                zt1 = fourpi/(omega*dble(nkptnr))
                wplas = zt1*dabs(wplas)
                if (rank==0) then
                    write(fname, '("PLASMA_", 2I1, ".OUT")') a, b
                    write(*, '("  Plasma frequency written to ", a)') trim(adjustl(fname))
                    open(60, File=trim(fname), Action='WRITE', Form='FORMATTED')
                    write(60, '(G18.10, " : plasma frequency")') sqrt(wplas)
                    close(60)
                end if
                zt1 = zi/fourpi
                do iw = 1, wgrid
                    sigma(iw) = sigma(iw)+zt1*wplas/(w(iw)+eta(iw))
                end do
            end if
        end if        

!-----------------------------------------------
! Print the output spectra
!-----------------------------------------------
        if (rank==0) then
! output energy units
            t1 = 1.0d0
            if (input%properties%dielmat%tevout) t1 = h2ev
! write the optical conductivity
            write(fname, '("SIGMA_", 2I1, ".OUT")') a, b
            write(*, '("  Optical conductivity tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            do iw = 1, wgrid
                write(60, '(3G18.10)') t1*w(iw), sigma(iw)
            end do
           close(60)
! write the dielectric function to file
            write(fname, '("EPSILON_", 2I1, ".OUT")') a, b
            write(*, '("  Dielectric tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            zt1 = zi*fourpi
            t2 = 0.0d0; if (a==b) t2 = 1.0d0
            do iw = 1, wgrid
                zt2 = t2+zt1*sigma(iw)/(w(iw)+eta(iw))
                write(60, '(3G18.10)') t1*w(iw), zt2
            end do
            close(60)
! write the EELS spectra
            write(fname, '("LOSS_", 2I1, ".OUT")') a, b
            write(*, '("  Electron energy loss spectrum written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            zt1 = zi*fourpi
            t2 = 0.0d0; if (a==b) t2 = 1.0d0
            do iw = 1, wgrid
                zt2 = t2+zt1*sigma(iw)/(w(iw)+eta(iw))
                write(60, '(2G18.10)') t1*w(iw), -aimag(1.0d0/zt2)
            end do
            close(60)
            if (input%properties%dielmat%tevout) write(*, '("  Output energy is in eV")')
            write(*,*)
        end if
#ifdef MPI        
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
    
    end do ! l, optical components
   
    close(50) 
    deallocate(pmat)
    deallocate(evalsvt,occsvt)
    deallocate(w,sigma,eta)
    
end subroutine
