!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine mygendmat ( lmax, fst, lst, is, ia, ngp1, ngp2, apwalm1, apwalm2, evecfv1, evecfv2, dmat)
      Use modinput
      Use modmain
      Implicit None
! arguments
      integer, intent( in) :: lmax, fst, lst, is, ia, ngp1, ngp2
      complex(8), intent( in) :: apwalm1( ngkmax, apwordmax, lmmaxapw, natmtot)
      complex(8), intent( in) :: apwalm2( ngkmax, apwordmax, lmmaxapw, natmtot)
      complex(8), intent( in) :: evecfv1( nmatmax, nstfv)
      complex(8), intent( in) :: evecfv2( nmatmax, nstfv)
      complex(8), intent( out) :: dmat( (lmax+1)**2, fst:lst, fst:lst)
! local variables
      integer :: lmmax
      integer :: l, m, lm
      integer :: i, j, n, ist, irc
      real(8) :: t1, t2
      complex(8) :: zt1
! automatic arrays
      real(8) :: fr1( nrcmtmax), fr2( nrcmtmax)
      real(8) :: gr( nrcmtmax), cf( 3, nrcmtmax)
! allocatable arrays
      complex(8), allocatable :: wfmt1(:,:), wfmt2(:,:)
      lmmax = (lmax+1)**2
! allocate local arrays
      allocate( wfmt1( lmmax, nrcmtmax))
      allocate( wfmt2( lmmax, nrcmtmax))
! zero the density matrix
      dmat = zzero
      n = lmmax * nrcmt (is)
! begin loop over second-variational states
      do j = fst, lst
        call wavefmt( input%groundstate%lradstep, lmax, is, ia, ngp1, apwalm1(:,:,:,:), evecfv1(:, j), lmmax, wfmt1)
        do i = fst, lst
          call wavefmt( input%groundstate%lradstep, lmax, is, ia, ngp2, apwalm2(:,:,:,:), evecfv2(:, i), lmmax, wfmt2)
          do l = 0, lmax
            do m = -l, l
              lm = idxlm( l, m)
              do irc = 1, nrcmt (is)
                zt1 = conjg( wfmt1( lm, irc))*wfmt2( lm, irc)
                t1 = rcmt( irc, is)**2
                fr1( irc) = dble( zt1)*t1
                fr2( irc) = aimag( zt1)*t1
              end do
              call fderiv( -1, nrcmt( is), rcmt( :, is), fr1, gr, cf)
              t1 = gr( nrcmt( is))
              call fderiv( -1, nrcmt( is), rcmt( :, is), fr2, gr, cf)
              t2 = gr( nrcmt( is))
              dmat( lm, j, i) = cmplx( t1, t2, 8)
            end do
          end do
        end do
      end do
      deallocate( wfmt1, wfmt2)
      return
End Subroutine
