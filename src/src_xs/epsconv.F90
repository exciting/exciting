


! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine epsconv
  use modmain
use modinput
  use modxs
  use modtetra
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='epsconv'
  character(256) :: filnam
  integer :: iq, iw, iwp, m, n, oct1, oct2, nc, un
  logical :: exis, tq0
  real(8), parameter :: epsc=1.d-8
  real(8), allocatable :: w(:), epst(:, :), lor(:), f(:), f1(:), g(:), g1(:), cf(:, :)
  complex(8), allocatable :: eps(:)
  integer, external :: l2int
  logical, external :: tqgamma
  ! initialize universal variables
  call init0
  call init1
  call init2
  ! original sampling method fo Brillouine zone
  bzsampl=l2int(input%xs%tetra%tetradf)
  allocate(w(nwdf), epst(nwdf, 2), eps(nwdf), lor(nwdf))
  allocate(f(nwdf), f1(nwdf), g(nwdf), g1(nwdf), cf(3, nwdf))
  ! loop over q-points
  do iq=1, nqpt
     ! matrix size for local field effects
     n=ngq(iq)
     tq0=tqgamma(iq)
     nc=1
     if (tq0) nc=3
     ! neglect/include local field effects
     do m=1, n, max(n-1, 1)
        ! loop over longitudinal components for optics
	do oct1=1, nc
	do oct2=1, nc
           ! generate filename for Tetrahedron method
	   call genfilname(basename = 'EPSILON', bzsampl = bzsampl, &
		nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),  fxctype = input%xs%tddft%fxctypenumber, &
		tq0 = tq0, oc1 = oct1, oc2 = oct2, iqmt = iq, filnam = filnam)
           ! check for file to read
	   inquire(file=trim(filnam), exist=exis)
	   if (.not.exis) then
	      write( * , * ) 'Error('//trim(thisnam)//'): file does not exist: '&
		   //trim(filnam)
	      call terminate
	   end if
	   call getunit(un)
	   open(unit = un, file = trim(filnam), form = 'formatted', &
		action = 'read', status = 'old')
           ! read energies and Re and Im of eps
	   do iw=1, nwdf
	      read(un, *) w(iw), epst(iw, 1), epst(iw, 2)
	   end do
	   close(un)
           ! generate filename for output (_s_ymmetric _l_orentzian)
	   call genfilname(basename = 'EPSILON_sl', bzsampl = bzsampl, &
		nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),  fxctype = input%xs%tddft%fxctypenumber, &
		tq0 = tq0, oc1 = oct1, oc2 = oct2, iqmt = iq, filnam = filnam)
	   open(unit = un, file = trim(filnam), &
		form = 'formatted', action = 'write', status = 'replace')
           ! rescale energies read from file to atomic units
	   w(:)=w(:)/escale
           ! complex dielectric function
	   eps(:)=epst(:, 1)+zi*epst(:, 2)
           ! convolution with Lorentzian
	   do iw=1, nwdf
	      do iwp=1, nwdf
!!$                 ! standard Lorentzian with peak at w(iw)
!!$                 lor(iwp)=(1/pi)*broad/((w(iw)-w(iwp))**2+broad**2)
                 ! antisymmetric Lorentzian at w(iw) and -w(iw)
                 ! with norm arctan(w/broad) to assure zero crossing
		 lor(iwp) = (1.d0/(2.d0 * atan(w(iw)/input%xs%broad))) * ( &
		      input%xs%broad/((w(iw) - w(iwp)) ** 2 + input%xs%broad** 2) - &
		      input%xs%broad/(( - w(iw) - w(iwp)) ** 2 + input%xs%broad** 2) )
		 if (w(iw) < epsc) lor(iwp) = 0.d0
		 f(iwp)=lor(iwp)*aimag(eps(iwp))
		 f1(iwp)=lor(iwp)*dble(eps(iwp))
	      end do
              ! do the convolution
	      call fderiv(-1, nwdf, w, f, g, cf)
	      call fderiv(-1, nwdf, w, f1, g1, cf)

	      write(un, '(4g18.10)') w(iw) * escale, g1(nwdf), g(nwdf), &
		   (pi/input%xs%broad) * input%xs%broad** 2/((w(iw) - w(iwp)) ** 2 + input%xs%broad** 2)
	   end do ! iw
!!$           call fsmooth(nsmdos,nwdf,1,g)
!!$           call fsmooth(nsmdos,nwdf,1,g1)
	   close(un)
	end do ! oct
	end do
     end do
  end do
  deallocate(w, epst, eps, lor, f, f1, g, g1, cf)
end subroutine epsconv
