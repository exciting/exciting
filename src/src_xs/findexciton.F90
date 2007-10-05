
module m_findoscstr
  implicit none
contains
  subroutine findoscstr(oct,nw,w,mdf1,ne,ei,o)
    use modmain
    use modtddft
    use modtetra
    implicit none
    ! arguments
    integer, intent(in) :: oct,nw,ne
    real(8), intent(in) :: w(:)
    integer, intent(in) :: ei(:)
    real(8), intent(in) :: mdf1(:)
    real(8) :: o(:)
    ! local variables
    real(8), allocatable :: f(:),fp(:),g(:),gp(:),cf(:,:),w2(:)
    integer :: iw,j,k,n
    
    allocate(f(nw),fp(nw),cf(3,nw),g(nw),gp(nw),w2(nw))
    f(:)=fxc0(1:nw,oct)
    g(:)=mdf1(1:nw)
    w2(:)=w(1:nw)

    ! derivative of RPA dielectric function
    call fderiv(1,nw,w2,g,gp,cf)
    ! derivative of xc-kernel
    call fderiv(1,nw,w2,f,fp,cf)

    do j=1,nexcit(oct)
       k=ei(j)
       o(j)=-pi/(f(k)*abs(fp(k)*(g(k)-1.d0)+f(k)*&
            gp(k)))
    end do

    deallocate(f,fp,cf,g,gp,w2)

  end subroutine findoscstr
end module m_findoscstr


module m_findroots
  implicit none
contains
  subroutine findroots(n,w,f,nr,r,ri)
    implicit none
    ! arguments
    integer, intent(in) :: n
    real(8), intent(in) :: w(:),f(:)
    integer, intent(out) :: nr, ri(:)
    real(8), intent(out) :: r(:)
    ! local variables
    real(8) :: f0,f1
    integer :: i,j,k
    i=0
    do j=1,n-1
       f0=f(j)
       f1=f(j+1)
       ! intersection with zero happens
       if (f0*f1<0.d0) then
          i=i+1
          ! tendency of root to either of both grid points
          k=0
          if (0.5d0*(f0+f1)>0.d0) k=1
          r(i)=w(j+k)
          ri(i)=j+k
       end if
    end do
    nr=i
  end subroutine findroots
end module m_findroots



module m_findexciton
  implicit none
contains
  subroutine findexciton(oct,nw,w)
    use modmain
    use modtddft
    use modtetra
    use m_findroots
    use m_findoscstr
    implicit none
    ! arguments
    integer, intent(in) :: oct,nw
    real(8), intent(in) :: w(:)
    ! local variables
    character(*), parameter :: thisnam='findexciton'
    real(8), parameter :: dogap=1.d-4
    real(8), allocatable :: f(:)
    real(8) :: t1,wgap
    integer, allocatable :: exciti(:)
    integer :: iw,j,n,nwg

    ! find optical gap
    t1=maxval(mdfrpa(:,oct,2))*dogap
    do iw=1,nwdos
       if (abs(dble(mdfrpa(iw,oct,2))) > t1) exit
    end do
    nwg=iw
    wgap=w(nwg)
    write(unitout,'(a,f12.6,a,f12.6,a)') 'Info('//thisnam//'): optical &
         &bandgap from tetrahedron method at:',w(nwg),'  (',w(nwg)*h2ev,'eV)'

    allocate(f(nwg),exciti(nwg))
    exciti(:)=0

    ! find roots of constituing equation for exciton energies
    f(:)=1.d0+fxc0(:,oct)*(mdfrpa(:,oct,1)-1.d0)

    call findroots(nwg,w,f,nexcit(oct),excite(:,oct),exciti)
    call findoscstr(oct,nwg,w,mdfrpa(:,oct,1),nexcit(oct),exciti,&
         excito(1:nexcit(oct),oct))

    if (nexcit(oct)>0) then
       write(unitout,'(a,2i6)') 'bound excitions found'
       write(unitout,'(a,2i6)') 'polarization, number of excitons:',oct,&
            nexcit(oct)
       write(unitout,'(a)') 'exciton energy (eV) and oscillator strength:'
       write(unitout,'(2f12.6)') (excite(j,oct)*h2ev,excito(j,oct),&
            j=1,nexcit(oct))
       write(unitout,'(a)') 'exciton binding energy (eV)'
       write(unitout,'(2f12.6)') ((wgap-excite(j,oct))*h2ev,&
            j=1,nexcit(oct))
    end if

    deallocate(f,exciti)

  end subroutine findexciton
end module m_findexciton

