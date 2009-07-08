

! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine getevecfv(vpl, vgpl, evecfv)
  use modmain
use modinput
  use modmpi
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  real(8), intent(in) :: vgpl(3, ngkmax)
  complex(8), intent(out) :: evecfv(nmatmax, nstfv, nspnfv)
  ! local variables
  logical::exist
  integer::isym, lspl, ilo, l, m, lm, koffset
  integer::ik, igp, igk, ig, i, ilspl
  integer::is, ia, ja, ias, jas
  integer::recl, nmatmax_, nstfv_, nspnfv_
  real(8)::vkl_(3), v(3), t1
  real(8)::si(3, 3), sc(3, 3)
  complex(8) zt1
  ! allocatable arrays
#ifdef XS
  ! added feature to access arrays for only a subset of bands
  complex(8), allocatable :: evecfv_(:, :, :)
#endif
  complex(8), allocatable :: evecfvt(:, :)
  complex(8), allocatable :: zflm1(:, :), zflm2(:, :)
  character(256) ::filetag
  character(256), external:: outfilenamestring
  ! find the equivalent k-point number and crystal symmetry element
  call findkpt(vpl, isym, ik)
  ! index to spatial rotation in lattice point group
  lspl=lsplsymc(isym)
  ! find the record length
#ifdef XS
  inquire(iolength=recl) vkl_, nmatmax_, nstfv_, nspnfv_
#endif
#ifndef XS
  inquire(iolength=recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
#endif
  !$OMP CRITICAL
  filetag=trim(filetag_evecfv)
  do i=1, 100
inquire(file=outfilenamestring(filetag, ik), exist=exist)
 if (exist) then
open(70, file = outfilenamestring(filetag, ik), action = 'READ', &
 form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
exit
else
call system('sync')
call sleep(5)
write(*, *) "Waiting for other process to write"//":getevecfv:"// &
     trim(outfilenamestring(filetag, ik))
endif
enddo
 if (splittfile) then
 koffset=ik-firstk(procofk(ik))+1
 else
 koffset =ik
 endif
#ifdef XS
  read(70, rec=1) vkl_, nmatmax_, nstfv_, nspnfv_
  close(70)
  if (nstfv.gt.nstfv_) then
     write(*, *)
     write(*, '("Error(getevecfv): invalid nstfv for k-point ", I8)') ik
     write(*, '(" current    : ", I8)') nstfv
     write(*, '(" EVECFV.OUT : ", I8)') nstfv_
     write(*, '(" file	     : ", a	 )') trim(outfilenamestring(filetag, ik))
     write(*, *)
     stop
  end if
  allocate(evecfv_(nmatmax_, nstfv_, nspnfv_))
  inquire(iolength=recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv_
  open(70, file = outfilenamestring(filetag, ik), action = 'READ', &
       form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
  read(70, rec=koffset) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv_
  ! retreive subset
  evecfv(:, :, :)=evecfv_(:, :nstfv, :)
  deallocate(evecfv_)
#endif
#ifndef XS
read(70, rec=koffset) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
#endif
  close(70)
  !$OMP END CRITICAL
  t1=abs(vkl(1, ik)-vkl_(1))+abs(vkl(2, ik)-vkl_(2))+abs(vkl(3, ik)-vkl_(3))
  if (t1.gt.input%structure%epslat) then
  write(*, *)
  write(*, '("Error(getevecfv): differing vectors for k-point ", I8)') ik
  write(*, '(" current	  : ", 3G18.10)') vkl(:, ik)
  write(*, '(" EVECFV.OUT : ", 3G18.10)') vkl_
  write(*, '(" file	  : ", a      )') trim(outfilenamestring(filetag, ik))
  write(*, *)
  stop
  end if
  if (nmatmax.ne.nmatmax_) then
  write(*, *)
  write(*, '("Error(getevecfv): differing nmatmax for k-point ", I8)') ik
  write(*, '(" current	  : ", I8)') nmatmax
  write(*, '(" EVECFV.OUT : ", I8)') nmatmax_
  write(*, '(" file	  : ", a      )') trim(outfilenamestring(filetag, ik))
  write(*, *)
  stop
  end if
#ifndef XS
  if (nstfv.ne.nstfv_) then
  write(*, *)
  write(*, '("Error(getevecfv): differing nstfv for k-point ", I8)') ik
  write(*, '(" current	  : ", I8)') nstfv
  write(*, '(" EVECFV.OUT : ", I8)') nstfv_
  write(*, '(" file	  : ", a      )') trim(outfilenamestring(filetag, ik))
  write(*, *)
  stop
  end if
#endif
  if (nspnfv.ne.nspnfv_) then
  write(*, *)
  write(*, '("Error(getevecfv): differing nspnfv for k-point ", I8)') ik
  write(*, '(" current	  : ", I8)') nspnfv
  write(*, '(" EVECFV.OUT : ", I8)') nspnfv_
  write(*, '(" file	  : ", a      )') trim(outfilenamestring(filetag, ik))
  write(*, *)
  stop
  end if
  ! if p = k then return
  t1=abs(vpl(1)-vkl(1, ik))+abs(vpl(2)-vkl(2, ik))+abs(vpl(3)-vkl(3, ik))
  if (t1.lt.input%structure%epslat) return
  if (input%groundstate%spin%spinsprl) then
  write(*, *)
  write(*, '("Error(getevec): code limitation - cannot rotate spin-spiral &
  &states")')
  write(*, '(" (first run one self-consistent loop with no k-point reduction)")')
  write(*, *)
  stop
  end if
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:, :)=symlat(:, :, ilspl)
!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
allocate(evecfvt(nmatmax, nstfv))
do igk=1, ngk(1, ik)
  ig=igkig(igk, 1, ik)
  v(:)=dble(ivg(:, ig))
  t1=-twopi*dot_product(v(:), vtlsymc(:, isym))
  zt1=cmplx(cos(t1), sin(t1), 8)
  evecfvt(igk, :)=zt1*evecfv(igk, :, 1)
end do
do igk=1, ngk(1, ik)
   call r3mtv(si, vgkl(:, igk, 1, ik), v)
   do igp=1, ngk(ik, 1)
      t1=abs(v(1)-vgpl(1, igp))+abs(v(2)-vgpl(2, igp))+abs(v(3)-vgpl(3, igp))
      if (t1.lt.input%structure%epslat) then
	 evecfv(igp, :, 1)=evecfvt(igk, :)
	 goto 10
      end if
   end do
10 continue
end do
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
if (nlotot.le.0) goto 20
allocate(zflm1(lolmmax, nstfv), zflm2(lolmmax, nstfv))
! make a copy of the local-orbital coefficients
do i=ngk(1, ik)+1, nmat(1, ik)
  evecfvt(i, :)=evecfv(i, :, 1)
end do
! spatial rotation symmetry matrix in Cartesian coordinates
sc(:, :)=symlatc(:, :, lspl)
! rotate k-point by inverse symmetry matrix
call r3mtv(si, vkl(:, ik), v)
do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia, is)
    ! equivalent atom for this symmetry
    ja=ieqatom(ia, is, isym)
    jas=idxas(ja, is)
    ! phase factor from translation
    t1 =- twopi * dot_product(vkl(:, ik), input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord(:))
    zt1=cmplx(cos(t1), sin(t1), 8)
    t1=twopi*dot_product(v(:), input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
    zt1=zt1*cmplx(cos(t1), sin(t1), 8)
    ! rotate local orbitals
    do ilo=1, nlorb(is)
      l=lorbl(ilo, is)
      zflm1(:, :)=0.d0
      do m=-l, l
	lm=idxlm(l, m)
	i=ngk(1, ik)+idxlo(lm, ilo, jas)
	zflm1(lm, :)=evecfvt(i, :)
      end do
      call rotzflm(sc, l, nstfv, lolmmax, zflm1, zflm2)
      do m=-l, l
	lm=idxlm(l, m)
	i=ngk(1, ik)+idxlo(lm, ilo, ias)
	evecfv(i, :, 1)=zt1*zflm2(lm, :)
      end do
    end do
  end do
end do
deallocate(zflm1, zflm2)
20 continue
deallocate(evecfvt)
return
end subroutine getevecfv

module m_getevecfvr
  implicit none
contains


subroutine getevecfvr(fname, isti, istf, vpl, vgpl, evecfv)
    use modmain
    implicit none
    ! arguments
    character(*), intent(in) :: fname
    integer, intent(in) :: isti, istf
    real(8), intent(in) :: vpl(3), vgpl(3, ngkmax)
    complex(8), intent(out) :: evecfv(:, :, :)
    ! local variables
    integer :: err
    complex(8), allocatable :: evecfvt(:, :, :)
    character(256) :: str
    ! check correct shapes
    err=0
    if ((isti.lt.1).or.(istf.gt.nstfv).or.(istf.le.isti)) then
       write(*, *)
       write(*, '("Error(getevecfvr): inconsistent limits for bands:")')
       write(*, '(" band limits  : ", 2i6)') isti, istf
       write(*, '(" maximum value: ", i6)') nstfv
       write(*, *)
       err=err+1
    end if
    if (size(evecfv, 2).ne.(istf-isti+1)) then
       write(*, *)
       write(*, '("Error(getevecfvr): output array does not match for bands:")')
       write(*, '(" band limits 	     : ", 2i6)') isti, istf
       write(*, '(" requested number of bands: ", i6)') istf-isti+1
       write(*, '(" array size		     : ", i6)') size(evecfv, 2)
       write(*, *)
       err=err+1
    end if
    if (size(evecfv, 1).ne.nmatmax) then
       write(*, *)
       write(*, '("Error(getevecfvr): output array does not match for &
	    &nmatmax:")')
       write(*, '(" nmatmax   : ", i6)') nmatmax
       write(*, '(" array size: ", i6)') size(evecfv, 1)
       write(*, *)
       err=err+1
    end if
    if (size(evecfv, 3).ne.nspnfv) then
       write(*, *)
       write(*, '("Error(getevecfvr): output array does not match for &
	    &nspnfv:")')
       write(*, '(" nspnfv    : ", i6)') nspnfv
       write(*, '(" array size: ", i6)') size(evecfv, 3)
       write(*, *)
       err=err+1
    end if
    if (err.ne.0) stop
    allocate(evecfvt(nmatmax, nstfv, nspnfv))
    filetag_evecfv=trim(fname)
    str=trim(filext)
    filext=''
    call getevecfv(vpl, vgpl, evecfvt)
    filetag_evecfv='EVECFV'
    filext=trim(str)
    evecfv(:, :, :)=evecfvt(:, isti:istf, :)
    deallocate(evecfvt)
  end subroutine getevecfvr
end module m_getevecfvr
