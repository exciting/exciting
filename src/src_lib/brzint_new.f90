!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: brzint_new
! !INTERFACE:
!
!
subroutine brzint_new( nsm, ngridk, nsk, ikmap, nw, wint, n, ld, e, f, g)
! !INPUT/OUTPUT PARAMETERS:
!   nsm    : level of smoothing for output function (in,integer)
!   ngridk : k-point grid size (in,integer(3))
!   nsk    : k-point subdivision grid size (in,integer(3))
!   ikmap  : map from grid to k-point set
!            (in,integer(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
!   nw     : number of energy divisions (in,integer)
!   wint   : energy interval (in,real(2))
!   n      : number of functions to integrate (in,integer)
!   ld     : leading dimension (in,integer)
!   e      : array of energies as a function of k-points (in,real(ld,*))
!   f      : array of weights as a function of k-points (in,real(ld,*))
!   g      : output function (out,real(nw))
! !DESCRIPTION:
!   See routine {\tt brzint}. This routine treats the integral over the trilinear
!   interpolated values analytically in two dimensions and numerically in the
!   third dimension (instead of numerically in all three dimensions).
!   This and some smaller changes makes this routine much faster and more accurate
!   for dense interpolation grids {\tt nsk} or a large number of energy divisions {\tt nw}.
!   The solution does usually not show strong noise.
!
! !REVISION HISTORY:
!   Created June 2018 (SeTi)
!EOP
!BOC
  implicit none
! arguments
  integer, intent( in) :: nsm, ngridk(3), nsk(3), ikmap(0:ngridk(1)-1, 0:ngridk(2)-1, 0:ngridk(3)-1), nw, n, ld
  real(8), intent( in) :: wint(2), e(ld, *), f( ld, *)
  real(8), intent( out) :: g( nw)
! local variables
  integer :: mat(8,8), i, j, ist, i3, j1, j2, j3, k1, k2, k3, iw, iw1, iw2, iint, inside(2,2), corner(2,2)
  real(8) :: z, wd, dw, dwi, w(nw), t1, ce(4), cf(4), a0, x0, x1, y0, y1, a, b, eps, intbound(2,2)
  real(8) :: e0(n,8), f0(n,8), ae(n,8), af(n,8)

  eps = 1.d-15

  mat(1,:) = (/ 1,  0,  0,  0,  0,  0,  0,  0/)
  mat(2,:) = (/-1,  1,  0,  0,  0,  0,  0,  0/)
  mat(3,:) = (/-1,  0,  1,  0,  0,  0,  0,  0/)
  mat(4,:) = (/-1,  0,  0,  0,  1,  0,  0,  0/)
  mat(5,:) = (/ 1, -1, -1,  1,  0,  0,  0,  0/)
  mat(6,:) = (/ 1, -1,  0,  0, -1,  1,  0,  0/)
  mat(7,:) = (/ 1,  0, -1,  0, -1,  0,  1,  0/)
  mat(8,:) = (/-1,  1,  1, -1,  1, -1, -1,  1/)

  if( minval( ngridk) .lt. 1) then
    write(*,*)
    write(*, '("Error(brzint): ngridk < 1 : ", 3I8)') ngridk
    write(*,*)
    stop
  end if
  if( minval( nsk) .lt. 1) then
    write(*,*)
    write(*, '("Error(brzint): nsk < 1 : ", 3I8)') nsk
    write(*,*)
    stop
  end if

! length of interval
  wd = wint(2) - wint(1)
! energy step size
  dw = wd/dble( nw - 1)
  do iw = 1, nw
    w( iw) = wint(1) + dble( iw - 1)*dw
  end do
  dwi = 1.d0/dw
  g = 0.d0


#ifdef USEOMP
!$omp parallel default( shared) private( j1, j2, j3, k1, k2, k3, e0, f0, ae, af, ist, i3, z, ce, cf, a, b, iw1, iw2, iw, a0, inside, corner, t1, x0, x1, y0, y1, iint, intbound, i, j)
!$omp do collapse(3)
#endif
  do j1 = 0, ngridk(1) - 1
    do j2 = 0, ngridk(2) - 1
      do j3 = 0, ngridk(3) - 1
        k1 = mod( j1+1, ngridk(1))
        k2 = mod( j2+1, ngridk(2))
        k3 = mod( j3+1, ngridk(3))

        e0(:,1) = e( 1:n, ikmap( j1, j2, j3))
        e0(:,2) = e( 1:n, ikmap( k1, j2, j3))
        e0(:,3) = e( 1:n, ikmap( j1, k2, j3))
        e0(:,4) = e( 1:n, ikmap( k1, k2, j3))
        e0(:,5) = e( 1:n, ikmap( j1, j2, k3))
        e0(:,6) = e( 1:n, ikmap( k1, j2, k3))
        e0(:,7) = e( 1:n, ikmap( j1, k2, k3))
        e0(:,8) = e( 1:n, ikmap( k1, k2, k3))

        f0(:,1) = f( 1:n, ikmap( j1, j2, j3))
        f0(:,2) = f( 1:n, ikmap( k1, j2, j3))
        f0(:,3) = f( 1:n, ikmap( j1, k2, j3))
        f0(:,4) = f( 1:n, ikmap( k1, k2, j3))
        f0(:,5) = f( 1:n, ikmap( j1, j2, k3))
        f0(:,6) = f( 1:n, ikmap( k1, j2, k3))
        f0(:,7) = f( 1:n, ikmap( j1, k2, k3))
        f0(:,8) = f( 1:n, ikmap( k1, k2, k3))

        call dgemm( 'n', 't', n, 8, 8, -1.d0, e0, n, dble( mat), 8, 0.d0, ae, n)
        call dgemm( 'n', 't', n, 8, 8,  1.d0, f0, n, dble( mat), 8, 0.d0, af, n)

        do ist = 1, n
          if( .not. (any( e0( ist, :) .lt. wint(1)-nw*wd) .or. (any( e0( ist, :) .gt. wint(2)+nw*wd)))) then
          
            do i3 = 0, nsk(3) - 1
              z = dble( i3)/dble( nsk(3))
              ce(1) = ae( ist, 1) + ae( ist, 4)*z
              ce(2) = ae( ist, 2) + ae( ist, 6)*z
              ce(3) = ae( ist, 3) + ae( ist, 7)*z
              ce(4) = ae( ist, 5) + ae( ist, 8)*z
              cf(1) = af( ist, 1) + af( ist, 4)*z
              cf(2) = af( ist, 2) + af( ist, 6)*z
              cf(3) = af( ist, 3) + af( ist, 7)*z
              cf(4) = af( ist, 5) + af( ist, 8)*z
              !write(*,'(5i,8g24.12e3)') j1, j2, j3, ist, i3, e0( ist, :)

              a = -maxval( (/ce(1), ce(1)+ce(2), ce(1)+ce(3), ce(1)+ce(2)+ce(3)+ce(4)/))
              b = -minval( (/ce(1), ce(1)+ce(2), ce(1)+ce(3), ce(1)+ce(2)+ce(3)+ce(4)/))
              iw1 = max( 1, ceiling( (a-wint(1))*dwi + 1.d0))
              iw2 = min( nw, floor( (b-wint(1))*dwi + 1.d0))

              if( (iw1 .ge. 1) .and. (iw2 .le. nw)) then

                do iw = iw1, iw2
                  a0 = ce(1) + w( iw)
                  inside = 0
                  corner = 0

                  if( abs( a0*ce(4) - ce(2)*ce(3)) .gt. eps) then
                    ! a3 != 0
                    if( abs( ce(4)) .gt. eps) then
                      t1 = 1.d0/ce(4)

                      x0 = -a0
                      if( abs( ce(2)) .gt. eps) then
                        x0 = x0/ce(2)
                      else
                        x0 = 2.d0
                      end if
                      if( (x0 .ge. eps) .and. (x0 .le. 1.d0-eps)) inside(1,1) = 1

                      x1 = -(a0 + ce(3))
                      if( abs( ce(2) + ce(4)) .gt. eps) then
                        x1 = x1/(ce(2) + ce(4))
                      else
                        x1 = 2.d0
                      end if
                      if( (x1 .ge. eps) .and. (x1 .le. 1.d0-eps)) inside(2,1) = 1
                      
                      y0 = -a0
                      if( abs( ce(3)) .gt. eps) then
                        y0 = y0/ce(3)
                      else
                        y0 = 2.d0
                      end if
                      if( (y0 .ge. eps) .and. (y0 .le. 1.d0-eps)) inside(1,2) = 1
                      if( (abs( y0) .lt. eps) .and. (abs( x0) .lt. eps)) corner(1,1) = 1
                      if( (abs( y0-1.d0) .lt. eps) .and. (abs( x1) .lt. eps)) corner(1,2) = 1
                      
                      y1 = -(a0 + ce(2))
                      if( abs( ce(3) + ce(4)) .gt. eps) then
                        y1 = y1/(ce(3) + ce(4))
                      else
                        y1 = 2.d0
                      end if
                      if( (y1 .ge. eps) .and. (y1 .le. 1.d0-eps)) inside(2,2) = 1
                      if( (abs( y1) .lt. eps) .and. (abs( x0-1.d0) .lt. eps)) corner(2,1) = 1
                      if( (abs( y1-1.d0) .lt. eps) .and. (abs( x1-1.d0) .lt. eps)) corner(2,2) = 1

                      iint = 0
                      ! inside cases
                      if( inside(1,2)) then
                        if( inside(2,2)) then
                          ! y0,y1,x0,x1 in
                          if( inside(1,1) .and. inside(2,1)) then
                            iint = 2
                            intbound(:,1) = (/0.d0, min( x0, x1)/)
                            intbound(:,2) = (/max( x0, x1), 1.d0/)
                          ! y0,y1 in, x0,x1 out
                          else if( .not. inside(1,1) .and. .not. inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, 1.d0/)
                          else
                            !write(*,*) "this should not happen 1"
                          end if
                        else
                          ! y0 in, y1 out, x0 out, x1 in
                          if( inside(2,1) .and. .not. inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x1/)
                          ! y0 in, y1 out, x0 in, x1 out
                          else if( inside(1,1) .and. .not. inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x0/)
                          else
                            !write(*,*) "this should not happen 2"
                          end if
                        end if
                      else
                        if( inside(2,2)) then
                          ! y0 out, y1 in, x0 in, x1 out
                          if( inside(1,1) .and. .not. inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/x0, 1.d0/)
                          ! y0 out, y1 in, x0 out, x1 in
                          else if( inside(2,1) .and. .not. inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/x1, 1.d0/)
                          else
                            !write(*,*) "this should not happen 3"
                          end if
                        else
                          ! y0,y1 out, x0,x1 in
                          if( inside(1,1) .and. inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/min( x0, x1), max( x0, x1)/)
                          !else
                          !  write(*,*) "this should not happen 4"
                          end if
                        end if
                      end if

                      ! single corner cases
                      if( corner(1,1)) then
                      !  write(*,'(5g24.12e3)') y0, x0, a0, ce(1,2), ce(2,1)
                        if( inside(2,2)) then
                          ! y0=0, y1 in, x1 in
                          if( inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/x1, 1.d0/)
                          ! y0=0, y1 in, x1 out
                          else
                            iint = 1
                            intbound(:,1) = (/0.d0, 1.d0/)
                          end if
                        else
                          ! y0=0, y1 out, x1 in
                          if( inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x1/)
                          end if
                        end if
                      else if( corner(1,2)) then
                        if( inside(2,2)) then
                          ! y0=1, y1 in, x0 in
                          if( inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/x0, 1.d0/)
                          ! y0=1, y1 in, x0 out
                          else
                            iint = 1
                            intbound(:,1) = (/0.d0, 1.d0/)
                          end if
                        else
                          !y0=1, y1 out, x0 in
                          if( inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x0/)
                          end if
                        end if
                      else if( corner(2,1)) then
                        if( inside(1,2)) then
                          ! y1=0, y0 in, x1 in
                          if( inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x1/)
                          ! y1=0, y0 in, x1 out
                          else
                            iint = 1
                            intbound(:,1) = (/0.d0, 1.d0/)
                          end if
                        else
                          ! y1=0, y0 out, x1 in
                          if( inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/x1, 1.d0/)
                          end if
                        end if
                      else if( corner(2,2)) then
                        if( inside(1,2)) then
                          ! y1=1, y0 in, x0 in
                          if( inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x0/)
                          ! y1=1, y0 in, x0 out
                          else
                            iint = 1
                            intbound(:,1) = (/0.d0, 1.d0/)
                          end if
                        else
                          ! y1=1, y0 out, x0 in
                          if( inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/x0, 1.d0/)
                          end if
                        end if
                      end if

                      do i = 1, iint
                        if( abs( intbound(2,i) - intbound(1,i)) .gt. eps) then
                          j = 1
                          if( abs( intbound(2,i)) .gt. abs( intbound(1,i))) j = 2
                          a = 0.d0
                          b = ce(4)*ce(4)*cf(1) - ce(3)*ce(4)*cf(2) - ce(2)*ce(4)*cf(3) - a0*ce(4)*cf(4) + 2.d0*ce(2)*ce(3)*cf(4)
                          a = a + b*log( abs( (ce(3)+ce(4)*intbound(2,i))/(ce(3)+ce(4)*intbound(1,i)))) 
                          b = ce(4)*ce(4)*cf(2) - ce(2)*ce(4)*cf(4)
                          a = a + b*(intbound(2,i) - intbound(1,i))
                          b = a0*ce(4)*ce(4)*cf(3) - ce(2)*ce(3)*ce(4)*cf(3) - a0*ce(3)*ce(4)*cf(4) + ce(2)*ce(3)*ce(3)*cf(4)
                          a = a + b*(1.d0/(ce(3) + ce(4)*intbound(2,i)) - 1.d0/(ce(3) + ce(4)*intbound(1,i)))
#ifdef USE  OMP
!$omp atomic update
#endif
                          if( abs( a) .gt. eps) g( iw) = g( iw) + t1*sign( t1*t1, ce(3)+ce(4)*intbound(j,i))*a
#ifdef USE  OMP
!$omp end atomic
#endif
                        end if
                      end do
                    
                    ! a3 = 0, a2 != 0
                    else if( abs( ce(3)) .gt. eps) then
                      !write(*,*) "apply a3=0 correction"
                      t1 = 1.d0/ce(3)

                      x0 = -a0
                      if( abs( ce(2)) .gt. eps) then
                        x0 = x0/ce(2)
                      else
                        x0 = 2.d0
                      end if
                      if( (x0 .ge. eps) .and. (x0 .le. 1.d0-eps)) inside(1,1) = 1

                      x1 = -(a0 + ce(3))
                      if( abs( ce(2)) .gt. eps) then
                        x1 = x1/ce(2)
                      else
                        x1 = 2.d0
                      end if
                      if( (x1 .ge. eps) .and. (x1 .le. 1.d0-eps)) inside(2,1) = 1
                      
                      y0 = -a0/ce(3)
                      if( (y0 .ge. eps) .and. (y0 .le. 1.d0-eps)) inside(1,2) = 1
                      
                      y1 = -(a0 + ce(2))/ce(3)
                      if( (y1 .ge. eps) .and. (y1 .le. 1.d0-eps)) inside(2,2) = 1
                      
                      if( inside(1,2)) then
                        if( inside(2,2)) then
                          iint = 1
                          intbound(:,1) = (/0.d0, 1.d0/)
                        else
                          if( inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x0/)
                          else if( inside(1,2)) then
                            iint = 1
                            intbound(:,1) = (/0.d0, x1/)
                          end if
                        end if
                      else
                        if( inside(2,2)) then
                          if( inside(1,1)) then
                            iint = 1
                            intbound(:,1) = (/x0, 1.d0/)
                          else if( inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/x1, 1.d0/)
                          end if
                        else
                          if( inside(1,1) .and. inside(2,1)) then
                            iint = 1
                            intbound(:,1) = (/min( x0, x1), max( x0, x1)/)
                          end if
                        end if
                      end if

                      do i = 1, iint
                        if( abs( intbound(2,i) - intbound(1,i)) .gt. eps) then
                          a = 0.d0
                          b = ce(3)*cf(1) - a0*cf(3)
                          a = a + b*(intbound(2,i) - intbound(1,i))
                          b = ce(3)*cf(2) - ce(2)*cf(3) - a0*cf(4)
                          a = a + b/2.d0*(intbound(2,i)**2 - intbound(1,i)**2)
                          b = ce(2)*cf(4)
                          a = a + b/3.d0*(intbound(2,i)**3 - intbound(1,i)**3)
#ifdef USE  OMP
!$omp atomic update
#endif
                          if( abs( a) .gt. eps) g( iw) = g( iw) + t1*abs( t1)*a
#ifdef USE  OMP
!$omp end atomic
#endif
                        end if
                      end do
                    else
                      write(*,*) "a3=a2=0 correction missing"
                    end if

                  end if
                  
                end do

              end if

            end do

          end if
        end do

      end do
    end do
  end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
! normalise function
  t1 = dble( ngridk(1)*ngridk(2)*ngridk(3))*dble( nsk(3))
  t1 = 1.d0/t1
  g = t1*g
! smooth output function if required
  if( nsm .gt. 0) call fsmooth( nsm, nw, 1, g)

  return
end subroutine
!EOC
