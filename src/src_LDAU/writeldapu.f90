
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writeldapu
use modmain
implicit none
! local variables
integer is,ia,ias,ispn,jspn
integer l,m1,m2,lm1,lm2
open(50,file='DMATLU'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
       trim(spsymb(is)),ia
      write(50,'(" l = ",I2)') l
      do ispn=1,nspinor
        do jspn=1,nspinor
          write(50,'(" ispn = ",I1,", jspn = ",I1)') ispn,jspn
          do m1=-l,l
            lm1=idxlm(l,m1)
            do m2=-l,l
              lm2=idxlm(l,m2)
              write(50,'("  m1 = ",I2,", m2 = ",I2," : ",2G18.10)') m1,m2, &
               dmatlu(lm1,lm2,ispn,jspn,ias)
            end do
          end do
        end do
      end do
    end do
  end if
end do
close(50)
open(50,file='VMATLU'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
       trim(spsymb(is)),ia
      write(50,'(" l = ",I2)') l
      do ispn=1,nspinor
        do jspn=1,nspinor
          write(50,'(" ispn = ",I1,", jspn = ",I1)') ispn,jspn
          do m1=-l,l
            lm1=idxlm(l,m1)
            do m2=-l,l
              lm2=idxlm(l,m2)
              write(50,'("  m1 = ",I2,", m2 = ",I2," : ",2G18.10)') m1,m2, &
               vmatlu(lm1,lm2,ispn,jspn,ias)
            end do
          end do
        end do
      end do
    end do
  end if
end do
close(50)
if (ldapu.eq.3) then
  open(50,file='ALPHALU'//trim(filext),action='WRITE',form='FORMATTED')
  do is=1,nspecies
    l=llu(is)
    if (l.ge.0) then
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        write(50,*)
        write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
         trim(spsymb(is)),ia
        write(50,'(" l = ",I2)') l
        write(50,'(" alpha = ",G18.10)') alphalu(ias)
      end do
    end if
  end do
  close(50)
end if
return
end subroutine

