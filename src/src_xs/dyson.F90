
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dyson
  implicit none
contains
  
  subroutine dyson(iq,oct,iw,n,s0,k,s)
    use invert
    !
    ! solve Dyson's equation S = S0 + S0*K*S
    ! for S by inversion.
    !
    implicit none
    ! arguments
    integer, intent(in) :: iq,oct,iw
    integer, intent(in) :: n
    complex(8), intent(in) :: s0(:,:), k(:,:)
    complex(8), intent(out) :: s(:,:)
    ! local variables
    character(*), parameter :: thisnam='dyson'
    complex(8),parameter :: zone=(1.d0,0.d0),zzero=(0.d0,0.d0)
    complex(8), allocatable :: mt(:,:)
    integer :: shs0(2),shk(2),shs(2),nmin,nmax,j


    integer :: a,b


    ! check matrix sizes
    shs0=shape(s0)
    shk=shape(k)
    shs=shape(s)
    nmin=minval((/shs0,shk,shs/))
    nmax=maxval((/shs0,shk,shs/))
    if ((nmin.ne.nmax).or.(nmin.lt.n)) then
       write(*,'("Error(",a,"): inconsistent matrix sizes")') trim(thisnam)
       write(*,'("  n :",i9)') n
       write(*,'("  S0:",2i9)') shs0
       write(*,'("  K :",2i9)') shk
       write(*,'("  S :",2i9)') shs
       call terminate
    end if

    ! allocate
    allocate(mt(n,n))

    ! calculate matrix -S0*K
    call zgemm('n','n', n, n, n, -zone, s0, n, k, n, zzero, mt, n )

    ! calculate matrix T=[1 - S0*K]
    forall(j=1:n) mt(j,j)=mt(j,j)+1.d0

    ! invert matrix T
    call zinvert_lapack(mt,mt)


    if ((iq.eq.1).and.(iw.eq.1).and.(n.ne.1)) then
       do a=1,n
          do b=1,n
             write(200+oct,'(2i8,2g18.10)') a,b,mt(a,b)
          end do
       end do


    end if


    ! calculate matrix S=T^-1*S0
    call zgemm('n','n', n, n, n, zone, mt, n, s0, n, zzero, s, n )

    deallocate(mt)

  end subroutine dyson
  


end module m_dyson



