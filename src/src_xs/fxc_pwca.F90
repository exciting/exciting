


! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

module m_fxc_alda
  implicit none
contains


subroutine fxc_alda(iq, msiz, fxcg)
    use modmain,only:ivg,gc,fourpi,lmmaxvr,nrmtmax,natmtot,ngrtot,ngvec,ivgig
    use modinput
    use modxs,only:unitout,fxcmt,gqc,fxcir,igqig
    use modtimer2
    use m_ftfun,only:ftfun
    implicit none
    ! arguments
    integer, intent(in) :: iq, msiz
    complex(8), intent(out) :: fxcg(:, :)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_alda'
    complex(8), allocatable :: fxcft1(:)
    integer :: sh(2), ig, ngmax, igq1, igq2, iv1(3), iv(3)
    sh=shape(fxcg)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout, '(a, 2i9, a, i9, a)') 'Error('//trim(thisnam)//'): size of &
	    &fxc is to small (required)', sh, '(', msiz, ')'
       stop
    end if
    if (allocated(fxcmt)) deallocate(fxcmt)
    if (allocated(fxcir)) deallocate(fxcir)
    allocate(fxcmt(lmmaxvr, nrmtmax, natmtot))
    allocate(fxcir(ngrtot))
    ! calculate exchange-correlation kernel in real space
    call kernxc
    ! determine G-vector cutoff for 2*|G+q|_max
    do ngmax=1, ngvec
       if (2.d0*input%xs%gqmax.lt.gc(ngmax)) exit
    end do
    ! Fourier transform of muffin-tin and interstitial kernel
    allocate(fxcft1(ngmax))
    call ftfun(ngmax, input%xs%tddft%lmaxalda, .true., .true., fxcir, fxcmt, fxcft1)
    ! transform G''=G-G' to G and G'
    do igq1=1, msiz
       iv1(:)=ivg(:, igqig(igq1, iq))
       do igq2=1, msiz
	  iv(:)=iv1(:)- ivg(:, igqig(igq2, iq))
	  ig = ivgig(iv(1), iv(2), iv(3))
	  fxcg(igq1, igq2)=fxcft1(ig)
          ! renormalization to symmetrized quantity wrt. G-space
	  fxcg(igq1, igq2)=fxcg(igq1, igq2)*(gqc(igq1, iq)*gqc(igq2, iq))/fourpi
       end do
    end do
    ! deallocate
    deallocate(fxcmt, fxcir, fxcft1)
  end subroutine fxc_alda

end module m_fxc_alda
