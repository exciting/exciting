
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevecfv(vpl,vgpl,evecfv)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer isym,id(3),ilo,l,m,lm
integer ik,igp,igk,ist,i
integer is,ia1,ia2,ias1,ias2
integer recl,nmatmax_,nstfv_,nspnfv_
real(8) vkl_(3),v1(3),v2(3),t1
real(8) s(3,3),si(3,3),sc(3,3)
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: evecfvt(:,:)
complex(8), allocatable :: zflm(:,:)
! external functions
real(8) r3taxi,r3dot
external r3taxi,r3dot
! find the k-point number
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nmatmax_,nstfv_,nspnfv_,evecfv
!$OMP CRITICAL
open(70,file=trim(scrpath)//'EVECFV'//trim(filext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nmatmax_,nstfv_,nspnfv_,evecfv
close(70)
!$OMP END CRITICAL
if (r3taxi(vkl(1,ik),vkl_).gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECFV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nmatmax.ne.nmatmax_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nmatmax for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nmatmax
  write(*,'(" EVECFV.OUT : ",I8)') nmatmax_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVECFV.OUT : ",I8)') nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVECFV.OUT : ",I8)') nspnfv_
  write(*,*)
  stop
end if
! if symmetry element is the identity return
if (isym.eq.1) return
if (spinsprl) then
  write(*,*)
  write(*,'("Error(getevec): code limitation - cannot rotate spin-spiral &
  &states")')
  write(*,'(" (first run one self-consistent loop with no k-point reduction)")')
  write(*,*)
  stop
end if
! real symmetry matrix
s(:,:)=dble(symcrys(:,:,isym))
! the inverse of s rotates k into p
call r3minv(s,si)
! translate and rotate APW coefficients
allocate(evecfvt(nmatmax,nstfv))
do ist=1,nstfv
  do igk=1,ngk(ik,1)
    evecfvt(igk,ist)=evecfv(igk,ist,1)
  end do
end do
do igk=1,ngk(ik,1)
  call r3mtv(si,vgkl(1,igk,ik,1),v1)
  do igp=1,ngk(ik,1)
    if (r3taxi(v1,vgpl(1,igp)).lt.epslat) then
      evecfv(igp,:,1)=evecfvt(igk,:)
      goto 10
    end if
  end do
10 continue
end do
deallocate(evecfvt)
! return if there are no local-orbitals
if (nlotot.le.0) return
! translate and rotate local-orbital coefficients
allocate(zflm(lolmmax,nstfv))
! convert symmetry matrix to Cartesian coordinates
call r3mm(s,ainv,sc)
call r3mm(avec,sc,sc)
do is=1,nspecies
  do ia1=1,natoms(is)
    ias1=idxas(ia1,is)
! find equivalent atom for this symmetry
    call r3mv(si,atposl(1,ia1,is),v1)
    call r3frac(epslat,v1,id)
    do ia2=1,natoms(is)
      v2(:)=atposl(:,ia2,is)
      call r3frac(epslat,v2,id)
      if (r3taxi(v1,v2).lt.epslat) goto 30
    end do
    write(*,*)
    write(*,'("Error(getevecfv): cannot transform one atom into another")')
    write(*,'(" with crystal symmetry : ",I4)') isym
    write(*,'(" for species ",I4)') is
    write(*,'(" and atom ",I4)') ia1
    write(*,*)
    stop
30 continue
    ias2=idxas(ia2,is)
! phase factor from translation
    t1=-twopi*r3dot(vkl(1,ik),atposl(1,ia1,is))
    zt1=cmplx(cos(t1),sin(t1),8)
    call r3mtv(si,vkl(1,ik),v1)
    t1=twopi*r3dot(v1,atposl(1,ia2,is))
    zt1=zt1*cmplx(cos(t1),sin(t1),8)
! rotate local orbitals
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      zflm(:,:)=0.d0
      do ist=1,nstfv
        do m=-l,l
          lm=idxlm(l,m)
          i=ngk(ik,1)+idxlo(lm,ilo,ias1)
          zflm(lm,ist)=evecfv(i,ist,1)
        end do
      end do
      call rotzflm(sc,l,nstfv,lolmmax,zflm,zflm)
      do ist=1,nstfv
        do m=-l,l
          lm=idxlm(l,m)
          i=ngk(ik,1)+idxlo(lm,ilo,ias2)
          evecfv(i,ist,1)=zt1*zflm(lm,ist)
        end do
      end do
    end do
  end do
end do
deallocate(zflm)
return
end subroutine

