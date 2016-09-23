subroutine findxslorb
  
  use modspecies

  implicit none
  integer :: nn
  real(8) :: e, de, t0, t1, dl, t00, t10, dl0
  real(8) :: p0(nrmt), p1(nrmt), q0(nrmt), q1(nrmt)
  real(8) :: eshift
  
  character(10) :: fname
  
! Energy step
  de = 0.0001d0

! Energy shift
  eshift = 0.0d0
  
  do l = 0, lxsmax
  
!   visualize the logarithmic derivative D_l
    write(fname,'("dl_l=",i1,".dat")') l
    write(*,*) fname
    open(777,file=fname,action='write')

    e = esccut
    call rschroddme(0, l, 0, e, nrmt, r, vr, nn, p0, p1, q0, q1)
    t00 = p0(nrmt)
    t10 = p1(nrmt)
    if ( dabs(t00)>1.0d-6 ) then
      dl0 = r(nrmt)*t10/t00+(l+1)
    else
      dl0 = 100.d0*dsign(1.d0,t00)*dsign(1.d0,t10)
    end if
    
!   search loop
    do while (e <= exsmax)
      
      e = e + de
      call rschroddme(0, l, 0, e, nrmt, r, vr, nn, p0, p1, q0, q1)
      t0 = p0(nrmt)
      t1 = p1(nrmt)
      
      if ( dabs(t0) > 1.d-6 ) then
        dl = r(nrmt)*t1/t0+(l+1)
        write(777,*) e, dl
      end if
      
!     if D_l changing sign
      if  ( dl*dl0 .le. 0.d0 ) then
        if ( dabs(dl).le.1.d0 ) then
          write(*,*) 'l=',l, '    El=', e
          
          ! treat only unoccupied region
          if (e>exscut) then
            nl(l) = nl(l)+1
            el(l,nl(l)) = e+eshift
          end if
        
        end if
        dl0 = dl
      end if

    end do ! e
    close(777)
    
  end do

end subroutine
