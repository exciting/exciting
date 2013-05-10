!
subroutine moke
    use modmain
    use modinput
    use modmpi
    implicit none

! local variables
    integer :: iw, iostat
    real(8) :: t1
    complex(8) :: exx, exy, zt1
    integer :: wgrid

! allocatable arrays
    real(8), allocatable :: w(:)
    real(8), allocatable :: eps11(:,:)
    real(8), allocatable :: eps12(:,:)
    complex(8), allocatable :: kerr(:)

! initialization for DIELMAT
    if (.not.associated(input%properties%dielmat)) &
   &  input%properties%dielmat => getstructdielmat(emptynode)
    
    input%properties%dielmat%intraband = input%properties%moke%intraband
    wgrid = input%properties%moke%wgrid
    input%properties%dielmat%wgrid = wgrid
    input%properties%dielmat%wmax = input%properties%moke%wmax
    input%properties%dielmat%scissor = input%properties%moke%scissor
    input%properties%dielmat%swidth = input%properties%moke%swidth
    input%properties%dielmat%drude = input%properties%moke%drude
    input%properties%dielmat%tevout = .false.

! for MOKE only two optical components required
    deallocate(input%properties%dielmat%epscomp)
    allocate(input%properties%dielmat%epscomp(2,2))
!   XX    
    input%properties%dielmat%epscomp(1,1) = 1
    input%properties%dielmat%epscomp(2,1) = 1
!   XY
    input%properties%dielmat%epscomp(1,2) = 1
    input%properties%dielmat%epscomp(2,2) = 2
    
    call dielmat

    call barrier
    if (rank==0) then

! allocate local arrays
        allocate(w(wgrid))
        allocate(eps11(wgrid,2),eps12(wgrid,2))
        allocate(kerr(wgrid))

! read diagonal contribution to optical conductivity
        open(50, File='EPSILON_11.OUT', Action='READ', Status='OLD', Form='FORMATTED', IoStat=IoStat)
        if (iostat .ne. 0) then
            write (*,*)
            write (*, '("Error(moke): error opening EPSILON_11.OUT")')
            write (*,*)
            stop
        end if
        do iw = 1, wgrid
            read(50, *) w(iw), eps11(iw,1)
        end do
        read(50,*)
        do iw = 1, wgrid
            read (50, *) w(iw), eps11(iw,2)
        end do
        close(50)

! read off-diagonal contribution to optical conductivity
        open (50, File='EPSILON_12.OUT', Action='READ', Status='OLD', Form='FORMATTED', IoStat=IoStat)
        if (iostat .ne. 0) then
            write (*,*)
            write (*, '("Error(moke): error opening EPSILON_12.OUT")')
            write (*,*)
            stop
        end if
        do iw = 1, wgrid
            read(50, *) w(iw), eps12(iw,1)
        end do
        read(50,*)
        do iw = 1, wgrid
            read(50, *) w(iw), eps12(iw,2)
        end do
        close(50)

! calculate the complex Kerr angle
        do iw = 1, wgrid
            exx = cmplx(eps11(iw,1),eps11(iw,2),8)
            exy = cmplx(eps12(iw,1),eps12(iw,2),8)
            zt1 = (exx-zone)*sqrt(exx)
            if (abs(zt1) > 1.d-8) then
                kerr(iw) = -exy / zt1
            else
                kerr(iw) = zzero
            end if
        end do

! output energy units
        t1 = 1.0d0
        if (input%properties%moke%tevout) t1 = h2ev

! output results
        open(50, File='KERR.OUT', Action='WRITE', Form='FORMATTED')
        do iw = 1, wgrid
            write(50, '(2G18.10)') t1*w(iw), dble(kerr(iw))*180.d0/pi
        end do
        write(50, *)
        do iw = 1, wgrid
            write(50, '(2G18.10)') t1*w(iw), aimag(kerr(iw))*180.d0/pi
        end do
        close(50)
        write(*,*)
        write(*, '("Info(moke):")')
        write(*, '(" complex Kerr angle in degrees written to KERR.OUT")')
        if (input%properties%moke%tevout) &
       &  write(*, '(" Output energy is in eV")')
        write(*,*)
        
        deallocate(w, eps11, eps12, kerr)
    
    end if ! rank
    
    return

end subroutine
