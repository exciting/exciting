
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seitzgen(hall,ngen,srgen,stgen)
implicit none
character(20), intent(in) :: hall
integer, intent(out) :: ngen
real(8), intent(out) :: srgen(3,3,12)
real(8), intent(out) :: stgen(3,12)
! local variables
logical pr
integer i,m,n,no,nop
integer axis,id(3)
! zero vector tolerance
real(8), parameter :: eps=1.d-6
real(8) av(3),r(3,3),v1(3),v2(3),v3(3)
character(20) str1,str2,str3
! external functions
real(8) r3taxi
external r3taxi
str1=trim(adjustl(hall))//' '
no=0
nop=0
axis=0
n=0
10 continue
! check for origin shift vector
if (scan(str1,'(').eq.1) then
  if (index(str1,'(0 0 1)').ne.0) then
    v1(1)=0.d0; v1(2)=0.d0; v1(3)=1.d0
  else if (index(str1,'(0 0 -1)').ne.0) then
    v1(1)=0.d0; v1(2)=0.d0; v1(3)=-1.d0
  else
    write(*,*)
    write(*,'("Error(seitzgen): origin-shift not available : ",A)') trim(str1)
    write(*,*)
    stop
  end if
  v1(:)=v1(:)/12.d0
! apply vector shift to all Seitz matrices
  do i=1,ngen
    v3(:)=-v1(:)
    call r3mv(srgen(:,:,i),v3,v2)
    v2(:)=v2(:)+stgen(:,i)
    stgen(:,i)=v2(:)+v1(:)
  end do
  goto 20
end if
m=scan(str1,' ')
if (m.le.1) goto 20
str2=str1(1:m-1)
n=n+1
!------------------------------!
!     lattice translations     !
!------------------------------!
if (n.eq.1) then
  stgen(:,1)=0.d0
  if (scan(str2,'P').ne.0) then
    ngen=1
  else if (scan(str2,'A').ne.0) then
    stgen(1,2)=0.d0
    stgen(2,2)=0.5d0
    stgen(3,2)=0.5d0
    ngen=2
  else if (scan(str2,'B').ne.0) then
    stgen(1,2)=0.5d0
    stgen(2,2)=0.d0
    stgen(3,2)=0.5d0
    ngen=2
  else if (scan(str2,'C').ne.0) then
    stgen(1,2)=0.5d0
    stgen(2,2)=0.5d0
    stgen(3,2)=0.d0
    ngen=2
  else if (scan(str2,'I').ne.0) then
    stgen(:,2)=0.5d0
    ngen=2
  else if (scan(str2,'R').ne.0) then
    stgen(1,2)=0.6666666666666666667d0
    stgen(2,2)=0.3333333333333333333d0
    stgen(3,2)=0.3333333333333333333d0
    stgen(1,3)=0.3333333333333333333d0
    stgen(2,3)=0.6666666666666666667d0
    stgen(3,3)=0.6666666666666666667d0
    ngen=3
  else if (scan(str2,'S').ne.0) then
    stgen(1,2)=0.3333333333333333333d0
    stgen(2,2)=0.3333333333333333333d0
    stgen(3,2)=0.6666666666666666667d0
    stgen(1,3)=0.6666666666666666667d0
    stgen(2,3)=0.6666666666666666667d0
    stgen(3,3)=0.3333333333333333333d0
    ngen=3
  else if (scan(str2,'T').ne.0) then
    stgen(1,2)=0.3333333333333333333d0
    stgen(2,2)=0.6666666666666666667d0
    stgen(3,2)=0.3333333333333333333d0
    stgen(1,3)=0.6666666666666666667d0
    stgen(2,3)=0.3333333333333333333d0
    stgen(3,3)=0.6666666666666666667d0
    ngen=3
  else if (scan(str2,'F').ne.0) then
    stgen(1,2)=0.d0
    stgen(2,2)=0.5d0
    stgen(3,2)=0.5d0
    stgen(1,3)=0.5d0
    stgen(2,3)=0.d0
    stgen(3,3)=0.5d0
    stgen(1,4)=0.5d0
    stgen(2,4)=0.5d0
    stgen(3,4)=0.d0
    ngen=4
  else
    write(*,*)
    write(*,'("Error(seitzgen): Lattice symbol ''",A,"'' not found")') &
     trim(str2)
    write(*,*)
    stop
  end if
! set the rotations to the identity
  do i=1,ngen
    srgen(1,1,i)=1.d0; srgen(1,2,i)=0.d0; srgen(1,3,i)=0.d0
    srgen(2,1,i)=0.d0; srgen(2,2,i)=1.d0; srgen(2,3,i)=0.d0
    srgen(3,1,i)=0.d0; srgen(3,2,i)=0.d0; srgen(3,3,i)=1.d0
  end do
! check if lattice is centrosymmetric
  if (scan(str2,'-').ne.0) then
    do i=ngen+1,2*ngen
      srgen(:,:,i)=-srgen(:,:,i-ngen)
      stgen(:,i)=stgen(:,i-ngen)
    end do
    ngen=2*ngen
  end if
end if
!-------------------------------!
!     rotation-translations     !
!-------------------------------!
if (n.ge.2) then
! determine if rotation is proper or improper
  if (scan(str2,'-').eq.1) then
    pr=.false.
! remove the minus sign
    str3=str2(2:)
    str2=str3
  else
    pr=.true.
  end if
! determine the order of rotation
  if (scan(str2,'1').eq.1) then
    no=1
  else if (scan(str2,'2').eq.1) then
    no=2
  else if (scan(str2,'3').eq.1) then
    no=3
  else if (scan(str2,'4').eq.1) then
    no=4
  else if (scan(str2,'6').eq.1) then
    no=6
  else
    write(*,*)
    write(*,'("Error(seitzgen): invalid rotation order for Hall symbol ''",A,&
     &"''")') trim(hall)
    write(*,*)
    stop
  end if
! determine the axis of rotation
  if (scan(str2,'x').ne.0) then
! a axis
    axis=1
  else if (scan(str2,'y').ne.0) then
! b axis
    axis=2
  else if (scan(str2,'z').ne.0) then
! c axis
    axis=3
  else if (scan(str2,'"').ne.0) then
! a+b
    axis=5
  else if (scan(str2,'*').ne.0) then
! a+b+c axis
    axis=6
  else if (n.eq.2) then
! default first rotation is along c
    axis=3
  else if ((n.eq.3).and.(no.eq.2)) then
! default second rotation
    if ((nop.eq.2).or.(nop.eq.4)) then
! a axis
      axis=1
    else if ((nop.eq.3).or.(nop.eq.6)) then
! a-b axis
      axis=4
    else
      write(*,*)
      write(*,'("Error(seitzgen): malformed Hall symbol ''",A,"''")') trim(hall)
      write(*,'(" for default second rotation")')
      write(*,*)
      stop
    end if
  else if ((n.eq.4).and.(no.eq.3)) then
! third rotation around a+b+c axis
    axis=6
  else if (no.eq.1) then
! arbitrary axis for identity
    axis=1
  else
    write(*,*)
    write(*,'("Error(seitzgen): malformed Hall symbol ''",A,"''")') trim(hall)
    write(*,*)
    stop
  end if
! determine axis vector
  av(:)=0.d0
  if (axis.eq.1) then
! a axis
    av(1)=1.d0
  else if (axis.eq.2) then
! b axis
    av(2)=1.d0
  else if (axis.eq.3) then
! c axis
    av(3)=1.d0
  else if (axis.eq.4) then
! a-b axis
    av(1)=1.d0
    av(2)=-1.d0
  else if (axis.eq.5) then
! a+b axis
    av(1)=1.d0
    av(2)=1.d0
  else if (axis.eq.6) then
! a+b+c axis
    av(:)=1.d0
  end if
! compute the rotation part of the Seitz matrix
  if (axis.eq.1) then
! a axis
    if (no.eq.1) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    else if (no.eq.2) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)=-1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)=-1.d0
    else if (no.eq.3) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 0.d0; r(2,3)=-1.d0
      r(3,1)= 0.d0; r(3,2)= 1.d0; r(3,3)=-1.d0
    else if (no.eq.4) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 0.d0; r(2,3)=-1.d0
      r(3,1)= 0.d0; r(3,2)= 1.d0; r(3,3)= 0.d0
    else if (no.eq.6) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)=-1.d0
      r(3,1)= 0.d0; r(3,2)= 1.d0; r(3,3)= 0.d0
    end if
  else if (axis.eq.2) then
! b axis
    if (no.eq.1) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    else if (no.eq.2) then
      r(1,1)=-1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)=-1.d0
    else if (no.eq.3) then
      r(1,1)=-1.d0; r(1,2)= 0.d0; r(1,3)= 1.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)=-1.d0; r(3,2)= 0.d0; r(3,3)= 0.d0
    else if (no.eq.4) then
      r(1,1)= 0.d0; r(1,2)= 0.d0; r(1,3)= 1.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)=-1.d0; r(3,2)= 0.d0; r(3,3)= 0.d0
    else if (no.eq.6) then
      r(1,1)= 0.d0; r(1,2)= 0.d0; r(1,3)= 1.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)=-1.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    end if
  else if (axis.eq.3) then
! c axis
    if (no.eq.1) then
      r(1,1)= 1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)= 1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    else if (no.eq.2) then
      r(1,1)=-1.d0; r(1,2)= 0.d0; r(1,3)= 0.d0
      r(2,1)= 0.d0; r(2,2)=-1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    else if (no.eq.3) then
      r(1,1)= 0.d0; r(1,2)=-1.d0; r(1,3)= 0.d0
      r(2,1)= 1.d0; r(2,2)=-1.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    else if (no.eq.4) then
      r(1,1)= 0.d0; r(1,2)=-1.d0; r(1,3)= 0.d0
      r(2,1)= 1.d0; r(2,2)= 0.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    else if (no.eq.6) then
      r(1,1)= 1.d0; r(1,2)=-1.d0; r(1,3)= 0.d0
      r(2,1)= 1.d0; r(2,2)= 0.d0; r(2,3)= 0.d0
      r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)= 1.d0
    end if
  else if (axis.eq.4) then
! a-b axis
    r(1,1)= 0.d0; r(1,2)=-1.d0; r(1,3)= 0.d0
    r(2,1)=-1.d0; r(2,2)= 0.d0; r(2,3)= 0.d0
    r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)=-1.d0
  else if (axis.eq.5) then
! a+b axis
    r(1,1)= 0.d0; r(1,2)= 1.d0; r(1,3)= 0.d0
    r(2,1)= 1.d0; r(2,2)= 0.d0; r(2,3)= 0.d0
    r(3,1)= 0.d0; r(3,2)= 0.d0; r(3,3)=-1.d0
  else if (axis.eq.6) then
! a+b+c axis
    r(1,1)= 0.d0; r(1,2)= 0.d0; r(1,3)= 1.d0
    r(2,1)= 1.d0; r(2,2)= 0.d0; r(2,3)= 0.d0
    r(3,1)= 0.d0; r(3,2)= 1.d0; r(3,3)= 0.d0
  end if
! check if axis is invariant with respect to rotation
  call r3mv(r,av,v1)
  if (r3taxi(av,v1).gt.eps) then
    write(*,*)
    write(*,'("Error(seitzgen): axis not invariant with respect to rotation")')
    write(*,'(" for Hall symbol ''",A,"''")') trim(hall)
    write(*,*)
    stop
  end if
! apply inverse for improper rotation
  if (.not.pr) r(:,:)=-r(:,:)
! increment Seitz matrix count
  ngen=ngen+1
! store rotation in main array
  srgen(:,:,ngen)=r(:,:)
! remove rotation symbol
  str3=str2(2:)
  str2=str3
! determine translations
  stgen(:,ngen)=0.d0
  if (scan(str2,'a').ne.0) then
    stgen(1,ngen)=stgen(1,ngen)+0.5d0
  end if
  if (scan(str2,'b').ne.0) then
    stgen(2,ngen)=stgen(2,ngen)+0.5d0
  end if
  if (scan(str2,'c').ne.0) then
    stgen(3,ngen)=stgen(3,ngen)+0.5d0
  end if
  if (scan(str2,'n').ne.0) then
    stgen(:,ngen)=stgen(:,ngen)+0.5d0
  end if
  if (scan(str2,'u').ne.0) then
    stgen(1,ngen)=stgen(1,ngen)+0.25d0
  end if
  if (scan(str2,'v').ne.0) then
    stgen(2,ngen)=stgen(2,ngen)+0.25d0
  end if
  if (scan(str2,'w').ne.0) then
    stgen(3,ngen)=stgen(3,ngen)+0.25d0
  end if
  if (scan(str2,'d').ne.0) then
    stgen(:,ngen)=stgen(:,ngen)+0.25d0
  end if
  if (scan(str2,'1').ne.0) then
    if (no.eq.3) then
      stgen(:,ngen)=stgen(:,ngen)+0.3333333333333333333d0*av(:)
    else if (no.eq.4) then
      stgen(:,ngen)=stgen(:,ngen)+0.25d0*av(:)
    else if (no.eq.6) then
      stgen(:,ngen)=stgen(:,ngen)+0.1666666666666666667d0*av(:)
    end if
  else if (scan(str2,'2').ne.0) then
    if (no.eq.3) then
      stgen(:,ngen)=stgen(:,ngen)+0.6666666666666666667d0*av(:)
    else if (no.eq.6) then
      stgen(:,ngen)=stgen(:,ngen)+0.3333333333333333333d0*av(:)
    end if
  else if (scan(str2,'3').ne.0) then
    if (no.eq.4) then
      stgen(:,ngen)=stgen(:,ngen)+0.75d0*av(:)
    end if
  else if (scan(str2,'4').ne.0) then
    if (no.eq.6) then
      stgen(:,ngen)=stgen(:,ngen)+0.6666666666666666667d0*av(:)
    end if
  else if (scan(str2,'5').ne.0) then
    if (no.eq.6) then
      stgen(:,ngen)=stgen(:,ngen)+0.8333333333333333333d0*av(:)
    end if
  end if
end if
str3=adjustl(str1(m:))
str1=str3
nop=no
goto 10
20 continue
! map translations to [0,1)
do i=1,ngen
  call r3frac(eps,stgen(:,i),id)
end do
return
end subroutine

