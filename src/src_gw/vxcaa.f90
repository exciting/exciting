
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: vxcaa
! !INTERFACE:
subroutine vxcaa(tapp,is,ia,ngp,apwalm,v,h)
! !USES:
    use modinput
    use modmain
    use modgw
! !INPUT/OUTPUT PARAMETERS:
!   tapp   : .true. if the Hamiltonian is to be applied to the input vector,
!            .false. if the full matrix is to be calculated (in,logical)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   v      : input vector to which H is applied if tapp is .true., otherwise
!            not referenced (in,complex(nmatmax))
!   h      : H applied to v if tapp is .true., otherwise it is the XC
!            matrix in packed form (inout,complex(npmatmax))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created August 2006 (RGA) (from hmlaa)
!EOP
!BOC
    implicit none
! arguments
    logical, intent(in) :: tapp
    integer(4), intent(in) :: is
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: ngp
    complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in) :: v(nmatmax)
    complex(8), intent(inout) :: h(*)
! local variables
    integer(4) :: ias,io1,io2
    integer(4) :: l1,l2,l3,m1,m2,m3,lm1,lm2,lm3
    complex(8) :: zt1,zsum
! automatic arrays
    complex(8) :: zv(ngp)
! external functions
    real(8)    :: polynom
    complex(8) :: zdotc
    external polynom, zdotc

    ias=idxas(ia,is)
    do l1=0,input%groundstate%lmaxmat
      do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do io1=1,apword(l1,is)
          zv(:)=0.d0
          do l3=0,input%groundstate%lmaxapw
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
              if (lm1.ge.lm3) then
                do io2=1,apword(l3,is)
                  zsum=0.d0
                  do l2=0,input%groundstate%lmaxvr
                    if (mod(l1+l2+l3,2).eq.0) then
                      do m2=-l2,l2
                        lm2=idxlm(l2,m2)
                        if ((l2.eq.0).or.(l1.ge.l3)) then
                          zt1=gntyry(lm1,lm2,lm3)*vxcraa(io1,l1,io2,l3,lm2,ias)
                        else
                          zt1=gntyry(lm1,lm2,lm3)*vxcraa(io1,l3,io2,l1,lm2,ias)
                        end if
                        zsum=zsum+zt1
                      end do
                    end if
                  end do
                  if (lm1.eq.lm3) zsum=zsum*0.5d0
                  if (abs(dble(zsum))+abs(aimag(zsum)).gt.1.d-14) then
                    call zaxpy(ngp,zsum,apwalm(1,io2,lm3,ias),1,zv,1)
                  end if
                end do
              end if
            end do
          end do
          call zmatinp(tapp,ngp,zone,apwalm(1,io1,lm1,ias),zv,v,h)
        end do
      end do
    end do
    return
end subroutine
!EOC

