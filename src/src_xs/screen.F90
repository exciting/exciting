


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine screen
  use modxs
use modinput
  use m_genfilname
  implicit none
  ! local variables
  integer :: nwdft
  nwdft=nwdf
  nwdf=1
  input%xs%emattype=1
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  ! call dielectric function with only one frequency point
  call df
  ! alternative for checking only:
!!$  call screendm
  nwdf=nwdft
  write(unitout, '(a)') "Info(screen): Screening finished"
end subroutine screen


subroutine screendm
  ! *
use modinput
  ! March 2008
  ! poor man's implementation for checking * do not use for production use
  ! *

  use modmain
  use modxs
  use m_genfilname
  use m_getpemat
  use m_findgntn0
  use m_xsgauntgen
  implicit none
  integer :: iq, ik, ikq, n, ist, jst, ig, igp
  complex(8), allocatable :: scrn(:, :)
  real(8), allocatable :: scis12(:, :)

!!$character(256) :: fname
!!$integer :: nstmin,nstmax
!!$integer, allocatable :: nstk(:)
!!$real(8),allocatable :: wvkl(:,:),we(:,:)

  call init0
  ! initialise universal variables
  call init1
  call xssave0
  call init2
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax), input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  call findgntn0(max(input%xs%lmaxapwwf, lolmax), max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
  call readfermi
  ! all band combinations
  input%xs%emattype=0
  call ematbdcmbs(input%xs%emattype)

  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3, nstsv, nstsv))
  allocate(deou(nstsv, nstsv), docc12(nstsv, nstsv), scis12(nstsv, nstsv))

!!$!*************************************************************************
!!$call wien2k_energy_init(112233,fname='energy_DM',nlheader=2)
!!$call wien2k_energy_inquire(fname=fname,nstmin=nstmin,nstmax=nstmax)
!!$call wien2k_energy_printinfo
!!$
!!$allocate(wvkl(3,nkpt),we(nstmax,nkpt),nstk(nkpt))
!!$
!!$write(444,*) vkc/abs(bvec(1,1))
!!$call wien2k_energy_fetch(112233,wvkl,we,nstk,vkc/abs(bvec(1,1)),epslat)
!!$
!!$do ik=1,nkpt
!!$write(300,'(i6,3g18.10)') ik,wvkl(:,ik)
!!$write(301,'(2i6)') ik,nstk(ik)
!!$do n=1,nstmax
!!$write(302,'(2i6,g18.10)') ik,n,we(n,ik)
!!$end do
!!$end do
!!$
!!$deallocate(wvkl,we,nstk)
!!$stop 'stephan stopped here'
!!$!*************************************************************************

  ! loop over q-points
  do iq=1, nqpt
     call init1offs(qvkloff(1, iq))
     n=ngq(iq)
     if (allocated(xiou)) deallocate(xiou)
     allocate(xiou(nstsv, nstsv, n))
     allocate(scrn(n, n))
     scrn(:, :)=zzero
     call ematrad(iq)
     call ematqalloc

     ! loop over k-points
     do ik=1, nkpt
	ikq=ikmapikq(ik, iq)
	call getdevaldoccsv(iq, ik, ikq, istl1, istu1, istl2, istu2, deou, docc12, &
	     scis12)
	call ematqk1(iq, ik)
        ! get matrix elements (exp. expr. or momentum op.)
	call getpemat(iq, ik, 'PMAT_SCR.OUT', '', m12=xiou, p12=pmou)
        ! absorb optical matrix elements (consider only 11-polarization)
	if (iq.eq.1) xiou(:, :, 1)=pmou(1, :, :)
        ! summation for dielectric matrix
	do ist=1, nstsv
	  do jst=1, nstsv
!********
if (ist.lt.jst) then

	    if (abs(docc12(ist, jst)).gt.input%groundstate%epsocc) then
	      do ig=1, n
		do igp=1, n
    scrn(ig, igp) = scrn(ig, igp) - docc12(ist, jst)/deou(ist, jst) * xiou(ist, jst, ig)* &
			   conjg(xiou(ist, jst, igp)) * wkpt(ik)/omega * 2.d0 ! **** 

		    end do
		 end do
	      end if
end if
!********
              ! end loop over band-combinations
	   end do
	end do
        ! end loop over k-points
     end do
     ! diagonal
     forall (ig=1:n) scrn(ig, ig)=scrn(ig, ig)+1.d0
     do ig=1, n
	do igp=1, n
	   write(1000+iq, '(2i6, 2g18.10)') ig, igp, scrn(ig, igp)
	end do
     end do    
     call ematqdealloc
     write(unitout, '(a, i8)') 'Info(screendm): dielectric matrix finished &
	  &for q - point:', iq
     deallocate(xiou, scrn)
     ! end loop over q-points
  end do
  call findgntn0_clear
  deallocate(deou, docc12, scis12, pmou)

end subroutine screendm
