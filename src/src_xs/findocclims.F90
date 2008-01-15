
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findocclims(iq,iocc0,iocc,iunocc0,iunocc,io0,io,iu0,iu)
  use modmain
  use modxs
  use m_genfilname
  ! arguments
  integer, intent(in) :: iq
  integer, intent(out) :: iocc0,iocc,iunocc0,iunocc
  integer, intent(out) :: io0(nkpt),io(nkpt),iu0(nkpt),iu(nkpt)
  ! local variables
  integer :: ik,ikq,i0,i
  do ik=1,nkpt
     ! k+q-point set
     ikq=ik
     if (iq.ne.0) ikq=ikmapikq(ik,iq)
     call genfilname(iq=iq,setfilext=.true.)

write(*,*) 'ik,ikq',ik,ikq
write(*,*) 'vkl(1,ikq)',vkl(:,ikq)
write(*,*) 'vkl0(1,ik)',vkl0(:,ik)
     call getoccsv(vkl(1,ikq),occsv(1,ikq))
     do i=1,nstsv
        if (occsv(i,ik).lt.epsocc) exit
     end do
     io(ik)=i-1
     do i=nstsv,1,-1
        if (occsv(i,ik).gt.(occmax-epsocc)) exit
     end do
     iu(ik)=i+1
     if (iq.ne.0) then
        ! k-point set (q=0)
write(*,*) 'ik,ikq',ik,ikq
write(*,*) 'vkl(1,ikq)',vkl(:,ikq)
write(*,*) 'vkl0(1,ik)',vkl0(:,ik)
        call genfilname(iq=0,setfilext=.true.)
        call getoccsv0(vkl0(1,ik),occsv0(1,ik))
        call genfilname(iq=iq,setfilext=.true.)
        do i0=1,nstsv
           if (occsv0(i0,ik).lt.epsocc) exit
        end do
        io0(ik)=i0-1
        do i0=nstsv,1,-1
           if (occsv0(i0,ik).gt.(occmax-epsocc)) exit
        end do
        iu0(ik)=i0+1
     else
        io0(ik)=io(ik)
        iu0(ik)=iu(ik)
     end if
  end do
  ! overall highest (partially) occupied state
  iocc0=maxval(io0)
  iocc=maxval(io)
  ! overall lowest (partially) unoccupied state
  iunocc0=minval(iu0)
  iunocc=minval(iu)

  ! *** assign nstval and nstcon ***
  nstval=max(iocc0,iocc)
  nstcon=nstsv-nstval

  if ((iocc0.ge.iunocc).or.(iocc.ge.iunocc0)) then
     write(*,'(a)') 'Info(findocclims): partially occupied states present'
  end if

  ! debug output
  if (dbglev.gt.0) then
     write(*,'(a)') 'Debug(findocclims):'
     write(*,'(a)') ' iocc0,iocc,iunocc0,iunocc below:'
     write(*,'(4i8)') iocc0,iocc,iunocc0,iunocc
     write(*,'(a)') ' ik,io0,iu,diff,io,iu0,diff below:'
     do ik=1,nkpt
        write(*,'(7i8)') ik,io0(ik),iu(ik),iu(ik)-io0(ik), &
             io(ik),iu0(ik),iu0(ik)-io(ik)
     end do
     write(*,*)
  end if

end subroutine findocclims
