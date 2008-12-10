
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

module m_fxc_alda
  implicit none
contains

  subroutine fxc_alda(iq,msiz,fxcg)
    use modmain
    use modxs
    use modtimer2
    use m_ftfun
    implicit none
    ! arguments
    integer, intent(in) :: iq,msiz
    complex(8), intent(out) :: fxcg(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_alda'
    complex(8), allocatable :: fxcft1(:)
    integer :: sh(2),ig,igq1,igq2,iv1(3),iv(3)
    type(timer2) :: t
    call new_timer2(t)
    call tic_timer2(t)
    sh=shape(fxcg)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout,'(a,2i9,a,i9,a)') 'Error('//trim(thisnam)//'): size of &
            &fxc is to small (required)', sh, '(', msiz, ')'
       stop
    end if
    if (allocated(fxcmt)) deallocate(fxcmt)
    if (allocated(fxcir)) deallocate(fxcir)
    allocate(fxcmt(lmmaxvr,nrmtmax,natmtot))
    allocate(fxcir(ngrtot))
    call tic_timer2(t)
    call report_timer2(t,unitout,'initialization')
    ! calculate exchange-correlation kernel in real space
    call kernxc
    call tic_timer2(t)
    call report_timer2(t,unitout,'kernxc')
    ! Fourier transform of muffin-tin and interstitial kernel
    allocate(fxcft1(ngvec))
    call ftfun(ngvec,.true.,.true.,fxcir,fxcmt,fxcft1)
    call tic_timer2(t)
    call report_timer2(t,unitout,'Fourier transform')
    ! transform G''=G-G' to G and G'
    do igq1=1,msiz
       iv1(:)=ivg(:,igqig(igq1,iq))
       do igq2=1,msiz
          iv(:)=iv1(:)- ivg(:,igqig(igq2,iq))
          ig = ivgig(iv(1),iv(2),iv(3))
          fxcg(igq1,igq2)=fxcft1(ig)
          ! renormalization to symmetrized quantity wrt. G-space
          fxcg(igq1,igq2)=fxcg(igq1,igq2)*(gqc(igq1,iq)*gqc(igq2,iq))/fourpi
       end do
    end do
    ! deallocate
    deallocate(fxcmt,fxcir,fxcft1)
    call toc_timer2(t)
    call report_timer2(t,unitout,'settting up fxc_GGp')
  end subroutine fxc_alda

end module m_fxc_alda
