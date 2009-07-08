

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine pmatrad
  use modmain
use modinput
  use modxs
  implicit none
  ! local variables
  integer :: is, ia, ias, nr, ir
  integer :: l1, m1, lm1, l3, m3, lm3
  integer :: ilo, ilo1, ilo2, io, io1, io2, j
  ! automatic arrays
  real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3, nrmtmax)
  ! allocatable arrays
  real(8), allocatable :: fapw(:), flo(:), dapwfr(:, :, :, :, :), dlofr(:, :, :, :, :)
  ! allocate local arrays for radial derivatives
  allocate(fapw(nrmtmax))
  allocate(dapwfr(lmmaxapw, nrmtmax, 3, apwordmax, lmmaxapw))
  dapwfr(:, :, :, :, :)=0.d0
  if (nlotot.gt.0) then
     allocate(flo(nrmtmax))
     allocate(dlofr(lmmaxapw, nrmtmax, 3, nlomax, -lolmax:lolmax))
     dlofr(:, :, :, :, :)=0.d0
  end if
  ripaa(:, :, :, :, :, :)=0.d0
  if (nlotot.gt.0) then
     ripalo(:, :, :, :, :, :)=0.d0
     riploa(:, :, :, :, :, :)=0.d0
     riplolo(:, :, :, :, :, :)=0.d0
  end if
  ! begin loop over species
  do is=1, nspecies
     nr=nrmt(is)
     do ir=1, nr
        ! calculate r^2
	r2(ir)=spr(ir, is)**2
     end do
     ! begin loop over atoms
     do ia=1, natoms(is)
	ias=idxas(ia, is)
        !--------------------!
        !     derivatives    !
        !--------------------!
        ! APW functions
	do l1=0, input%groundstate%lmaxapw
	   do m1=-l1, l1
	      lm1=idxlm(l1, m1)
	      do io=1, apword(l1, is)
		 fapw(:)=apwfr(:, 1, io, l1, ias)
		 call gradzfmtr(input%groundstate%lmaxapw, nr, spr(1, is), l1, m1, lmmaxapw, &
		      nrmtmax, fapw, dapwfr(1, 1, 1, io, lm1))
	      end do
	   end do
	end do
	if (nlotot.gt.0) then
           ! local orbital functions
	   do ilo=1, nlorb(is)
	      l1=lorbl(ilo, is)
	      do m1=-l1, l1
		 lm1=idxlm(l1, m1)
		 flo(:)=lofr(:, 1, ilo, ias)
		 call gradzfmtr(input%groundstate%lmaxapw, nr, spr(1, is), l1, m1, lmmaxapw, &
		      nrmtmax, flo, dlofr(1, 1, 1, ilo, m1))
	      end do
	   end do
	end if
        !----------------!
        !     APW-APW    !
        !----------------!
	do l1=0, input%groundstate%lmaxapw
	   do m1=-l1, l1
	      lm1=idxlm(l1, m1)
	      do io1=1, apword(l1, is)
		 do l3=0, input%groundstate%lmaxapw
		    do m3=-l3, l3
		       lm3=idxlm(l3, m3)
		       do io2=1, apword(l3, is)
			  do j=1, 3
			     fr(:) = r2(1:nr) * apwfr(1:nr, 1, io1, l1, ias)* &
				  dapwfr(lm1, 1:nr, j, io2, lm3)
			     call fderiv(-1, nr, spr(1, is), fr, gr, cf)
			     ripaa(io1, lm1, io2, lm3, ias, j)=gr(nr)
			  end do
		       end do
		    end do
		 end do
	      end do
	   end do
	end do
	if (nlotot.gt.0) then
           !----------------------------!
           !     APW-local-orbital      !
           !----------------------------!
	   do l1=0, input%groundstate%lmaxapw
	      do m1=-l1, l1
		 lm1=idxlm(l1, m1)
		 do io1=1, apword(l1, is)
		    do ilo=1, nlorb(is)
		       l3=lorbl(ilo, is)
		       do m3=-l3, l3
			  lm3=idxlm(l3, m3)
			  do j=1, 3
			     fr(:) = r2(1:nr) * apwfr(1:nr, 1, io1, l1, ias)* &
				  dlofr(lm1, 1:nr, j, ilo, m3)
			     call fderiv(-1, nr, spr(1, is), fr, gr, cf)
			     ripalo(io1, lm1, ilo, m3, ias, j)=gr(nr)
			  end do
		       end do
		    end do
		 end do
	      end do
	   end do
           !----------------------------!
           !     local-orbital-APW      !
           !----------------------------!
	   do ilo=1, nlorb(is)
	      l1=lorbl(ilo, is)
	      do m1=-l1, l1
		 lm1=idxlm(l1, m1)
		 do l3=0, input%groundstate%lmaxapw
		    do m3=-l3, l3
		       lm3=idxlm(l3, m3)
		       do io2=1, apword(l3, is)
			  do j=1, 3
			     fr(:) = r2(1:nr) * lofr(:, 1, ilo, ias)* &
				  dapwfr(lm1, 1:nr, j, io2, lm3)
			     call fderiv(-1, nr, spr(1, is), fr, gr, cf)
			     riploa(ilo, m1, io2, lm3, ias, j)=gr(nr)
			  end do
		       end do
		    end do
		 end do
	      end do
	   end do
           !------------------------------------!
           !     local-orbital-local-orbital    !
           !------------------------------------!
	   do ilo1=1, nlorb(is)
	      l1=lorbl(ilo1, is)
	      do m1=-l1, l1
		 lm1=idxlm(l1, m1)
		 do ilo2=1, nlorb(is)
		    l3=lorbl(ilo2, is)
		    do m3=-l3, l3
		       lm3=idxlm(l3, m3)
		       do j=1, 3
			  fr(:) = r2(1:nr) * lofr(:, 1, ilo1, ias)* &
			       dlofr(lm1, 1:nr, j, ilo2, m3)
			  call fderiv(-1, nr, spr(1, is), fr, gr, cf)
			  riplolo(ilo1, m1, ilo2, m3, ias, j)=gr(nr)
		       end do
		    end do
		 end do
	      end do
	   end do
	end if
        ! end loops over atoms and species
     end do
  end do
  ! deallocate
  deallocate(fapw, dapwfr)
  if (nlotot.gt.0) deallocate(flo, dlofr)
end subroutine pmatrad
