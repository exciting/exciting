


! Copyright (C) 2002-2008 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine dielectric
use modmain
use modinput
implicit none
! local variables
integer::ik, jk, isym
integer::ist, jst, iw, i, j, l
integer::recl, nsk(3)
real(8)::eji, wd(2), wplas, t1, t2
real(8)::v1(3), v2(3), v3(3)
complex(8) zv(3), eta, zt1
character(256)::fname
! allocatable arrays
integer, allocatable :: lspl(:)
real(8), allocatable :: w(:)
real(8), allocatable :: f(:, :)
real(8), allocatable :: delta(:, :, :)
complex(8), allocatable :: pmat(:, :, :)
complex(8), allocatable :: sigma(:)
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
do ik=1, nkpt
! get the eigenvalues and occupancies from file
  call getevalsv(vkl(:, ik), evalsv(:, ik))
  call getoccsv(vkl(:, ik), occsv(:, ik))
end do
! allocate local arrays
allocate(lspl(nkptnr))
allocate(w(input%properties%dos%nwdos))
if (input%xs%tddft%intraband) allocate(f(nstsv, nkpt))
if (input%xs%usegdft) allocate(delta(nstsv, nstsv, nkpt))
allocate(pmat(3, nstsv, nstsv))
allocate(sigma(input%properties%dos%nwdos))
! compute generalised DFT correction
if (input%xs%usegdft) then
  call readstate
  call poteff
  call linengy
  call genapwfr
  call genlofr
  do ik=1, nkpt
    call gdft(ik, delta(:, :, ik))
  end do
end if
! energy interval should start from zero
wdos(1)=0.d0
! generate energy grid
t1=(wdos(2)-wdos(1))/dble(input%properties%dos%nwdos)
do iw=1, input%properties%dos%nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do
! find crystal symmetries which map non-reduced k-points to reduced equivalents
do ik=1, nkptnr
  call findkpt(vklnr(:, ik), isym, jk)
  lspl(ik)=lsplsymc(isym)
end do
! find the record length for momentum matrix element file
inquire(iolength=recl) pmat
open(50, file='PMAT.OUT', action='READ', form='UNFORMATTED', access='DIRECT', &
 recl=recl)
! i divided by the complex relaxation time
eta=cmplx(0.d0, input%groundstate%swidth)
! loop over dielectric tensor components
do l=1, noptcomp
  i=input%properties%linresponsetensor%optcomp(1, l)
  j=input%properties%linresponsetensor%optcomp(2, l)
  sigma(:)=0.d0
  if (input%xs%tddft%intraband) f(:, :)=0.d0
! loop over non-reduced k-points
  do ik=1, nkptnr
! equivalent reduced k-point
    jk=ikmap(ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
! read momentum matrix elements from direct-access file
    read(50, rec=jk) pmat
! valance states
    do ist=1, nstsv
      if (evalsv(ist, jk).lt.efermi) then
! conduction states
	do jst=1, nstsv
	  if (evalsv(jst, jk).gt.efermi) then
! rotate the matrix elements from the reduced to non-reduced k-point
! (note that the inverse operation is used)
	    v1(:)=dble(pmat(:, ist, jst))
	    call r3mv(symlatc(:, :, lspl(ik)), v1, v2)
	    v1(:)=aimag(pmat(:, ist, jst))
	    call r3mv(symlatc(:, :, lspl(ik)), v1, v3)
	    zv(:)=cmplx(v2(:), v3(:), 8)
	    zt1=occmax*zv(i)*conjg(zv(j))
	    eji=evalsv(jst, jk)-evalsv(ist, jk)+input%properties%bandstructure%scissor
	    if (input%xs%usegdft) eji=eji+delta(jst, ist, jk)
	    t1=1.d0/(eji+input%groundstate%swidth)
	    do iw=1, input%properties%dos%nwdos
	      sigma(iw) = sigma(iw) + t1 * (zt1/(w(iw) - eji + eta) &
	       +conjg(zt1)/(w(iw) + eji + eta))
	    end do
	  end if
	end do
      end if
    end do
  end do
  zt1=zi/(omega*dble(nkptnr))
  sigma(:)=zt1*sigma(:)
! intraband contribution
  if (input%xs%tddft%intraband) then
    if (i.eq.j) then
! compute plasma frequency
      do ik=1, nkpt
	do ist=1, nstsv
	  zt1=pmat(i, ist, ist)
	  f(ist, ik)=dble(zt1)**2+aimag(zt1**2)
	end do
      end do
      wd(1)=efermi-input%groundstate%swidth
      wd(2)=efermi+input%groundstate%swidth
! number of subdivisions used for interpolation
      nsk(:)=max(input%properties%dos%ngrdos/input%groundstate%ngridk(:), 1)
      call brzint(0, input%groundstate%ngridk, nsk, ikmap, 1, wd, nstsv, nstsv, evalsv, f, wplas)
      wplas=abs(wplas)*occmax*4.d0*pi/omega
      wplas=sqrt(wplas)
! write the plasma frequency to file
      write(fname, '("PLASMA_", 2I1, ".OUT")') i, j
      open(60, file=trim(fname), action='WRITE', form='FORMATTED')
      write(60, '(G18.10, " : plasma frequency")') wplas
      close(60)
! add the intraband contribution to sigma
      t1=wplas**2/fourpi
      do iw=1, input%properties%dos%nwdos
	sigma(iw)=sigma(iw)+t1/(input%groundstate%swidth-zi*w(iw))
      end do
    end if
  end if
! write the optical conductivity to file
  write(fname, '("SIGMA_", 2I1, ".OUT")') i, j
  open(60, file=trim(fname), action='WRITE', form='FORMATTED')
  do iw=1, input%properties%dos%nwdos
    write(60, '(2G18.10)') w(iw), dble(sigma(iw))
  end do
  write(60, '("     ")')
  do iw=1, input%properties%dos%nwdos
    write(60, '(2G18.10)') w(iw), aimag(sigma(iw))
  end do
  close(60)
! write the dielectric function to file
  write(fname, '("EPSILON_", 2I1, ".OUT")') i, j
  open(60, file=trim(fname), action='WRITE', form='FORMATTED')
  t1=0.d0
  if (i.eq.j) t1=1.d0
  do iw=1, input%properties%dos%nwdos
    if (w(iw).gt.1.d-8) then
      t2=t1-fourpi*aimag(sigma(iw)/(w(iw)+eta))
      write(60, '(2G18.10)') w(iw), t2
    end if
  end do
  write(60, '("     ")')
  do iw=1, input%properties%dos%nwdos
    if (w(iw).gt.1.d-8) then
      t2=fourpi*dble(sigma(iw)/(w(iw)+eta))
      write(60, '(2G18.10)') w(iw), t2
    end if
  end do
  close(60)
! end loop over tensor components
end do
write(*, *)
write(*, '("Info(dielectric):")')
write(*, '(" dielectric tensor written to EPSILON_ij.OUT")')
write(*, '(" optical conductivity written to SIGMA_ij.OUT")')
if (input%xs%tddft%intraband) then
  write(*, '(" plasma frequency written to PLASMA_ij.OUT")')
end if
write(*, '(" for components")')
do l=1, noptcomp
  write( * , '("  i = ", I1, ", j = ", I1)') input%properties%linresponsetensor%optcomp(1:2, l)
end do
write(*, *)
deallocate(lspl, w, pmat, sigma)
if (input%xs%tddft%intraband) deallocate(f)
if (input%xs%usegdft) deallocate(delta)
return
end subroutine
