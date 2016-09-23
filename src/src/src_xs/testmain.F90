!
! Copyright (C) 2004-2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Subroutine testmain
      Use modmain
      Use modxs
      Implicit None
      integer, parameter :: nd=1
      integer, parameter :: n(nd)=(/3/)
      integer, parameter :: nn=n(1)
      integer :: nt(nd), j
      complex(8) :: c(nn), cbw(nn), cfw(nn)

  ! test goes here
      write(*,*)
      write(*,'("Info(testmain): this is testmain.")')
      write(*,*)

      nt(:)=n(:)
      do j=1,nd
        call nfftifc(nt(j))
      end do
      if (any(nt.ne.n)) stop 'not a valid FFT grid'

  !-----------------------------------------------------------------------------

  ! initialize f=1
      c(:)=1.0d0

  ! backward transform to real-space
      cbw(:)=c(:)
      Call zfftifc (nd, n, 1, cbw)

  ! forward transform to G-space
      cfw(:)=c(:)
      Call zfftifc (nd, n, -1, cfw)

  ! write out
      open(500,file='FT_f=1',action='write',status='replace')
      write(500,*) '# untransformed (original data)'
      write(500,'(2g18.10)') (c(j),j=1,nn)
      close(500)
      open(501,file='FT_f=1_to_real-space',action='write',status='replace')
      write(501,*) '# backward transformed (to real-space)'
      write(501,'(2g18.10)') (cbw(j),j=1,nn)
      close(501)
      open(502,file='FT_f=1_to_G-space',action='write',status='replace')
      write(502,*) '# forward transformed (to G-space)'
      write(502,'(2g18.10)') (cfw(j),j=1,nn)
      close(502)

  !-----------------------------------------------------------------------------

  ! initialize f=delta_(j,1), only G=0 coefficient equal 1
      c(:)=0.0d0
      c(1)=1.d0

  ! backward transform to real-space
      cbw(:)=c(:)
      Call zfftifc (nd, n, 1, cbw)

  ! forward transform to G-space
      cfw(:)=c(:)
      Call zfftifc (nd, n, -1, cfw)

  ! write out
      open(500,file='FT_f=delta',action='write',status='replace')
      write(500,*) '# untransformed (original data)'
      write(500,'(2g18.10)') (c(j),j=1,nn)
      close(500)
      open(501,file='FT_f=delta_to_real-space',action='write',status='replace')
      write(501,*) '# backward transformed (to real-space)'
      write(501,'(2g18.10)') (cbw(j),j=1,nn)
      close(501)
      open(502,file='FT_f=delta_to_G-space',action='write',status='replace')
      write(502,*) '# forward transformed (to G-space)'
      write(502,'(2g18.10)') (cfw(j),j=1,nn)
      close(502)

  !-----------------------------------------------------------------------------

  ! initialize f (custom values)
      if (nn.lt.3) stop 'FFT array size < 3'
      c(:)=0.0d0
      c(1)=1.0d0
      c(2)=3.0d0
      c(3)=7.0d0

  ! backward transform to real-space
      cbw(:)=c(:)
      Call zfftifc (nd, n, 1, cbw)

  ! forward transform to G-space
      cfw(:)=c(:)
      Call zfftifc (nd, n, -1, cfw)

  ! write out
      open(500,file='FT_f=custom',action='write',status='replace')
      write(500,*) '# untransformed (original data)'
      write(500,'(2g18.10)') (c(j),j=1,nn)
      close(500)
      open(501,file='FT_f=custom_to_real-space',action='write',status='replace')
      write(501,*) '# backward transformed (to real-space)'
      write(501,'(2g18.10)') (cbw(j),j=1,nn)
      close(501)
      open(502,file='FT_f=custom_to_G-space',action='write',status='replace')
      write(502,*) '# forward transformed (to G-space)'
      write(502,'(2g18.10)') (cfw(j),j=1,nn)
      close(502)

End Subroutine
