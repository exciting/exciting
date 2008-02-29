
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findocclims(iq,iocc0,iocc,iunocc0,iunocc,io0,io,iu0,iu)
  use modmain
  use modxs
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  integer, intent(out) :: iocc0,iocc,iunocc0,iunocc
  integer, intent(out) :: io0(nkpt),io(nkpt),iu0(nkpt),iu(nkpt)
  ! local variables
  integer :: ik,ikq,i0,i
  logical :: t
  t=allocated(evalsv0)
  if (.not.t) allocate(evalsv0(nstsv,nkpt))
  do ik=1,nkpt
     ! k+q-point set
     ikq=ik
     if (iq.ne.0) ikq=ikmapikq(ik,iq)
     call getoccsv(vkl(1,ikq),occsv(1,ikq))
     call getevalsv(vkl(1,ikq),evalsv(1,ikq))
     do i=1,nstsv
        if (occsv(i,ikq).lt.epsocc) exit
     end do
     io(ik)=i-1
     do i=nstsv,1,-1
        if (occsv(i,ikq).gt.(occmax-epsocc)) exit
     end do
     iu(ik)=i+1
     if (iq.ne.0) then
        ! k-point set (q=0)
        call getoccsv0(vkl0(1,ik),occsv0(1,ik))
        call getevalsv0(vkl0(1,ik),evalsv0(1,ik))
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
  if (iq.ne.0) then
     ! lowest and highest valence energy
     evlmin=min(minval(evalsv(1,:)),minval(evalsv0(1,:)))
     evlmax=max(maxval(evalsv(nstsv,:)),maxval(evalsv0(nstsv,:)))
     ! lower and higher cutoff valence energy
     evlmincut=max(maxval(evalsv(1,:)),maxval(evalsv0(1,:)))
     evlmaxcut=min(minval(evalsv(nstsv,:)),minval(evalsv0(nstsv,:)))
  else
     ! lowest and highest valence energy
     evlmin=minval(evalsv(1,:))
     evlmax=maxval(evalsv(nstsv,:))
     ! lower and higher cutoff valence energy
     evlmincut=maxval(evalsv(1,:))
     evlmaxcut=minval(evalsv(nstsv,:))
  end if
  ! overall highest (partially) occupied state
  iocc0=maxval(io0)
  iocc=maxval(io)
  ! overall lowest (partially) unoccupied state
  iunocc0=minval(iu0)
  iunocc=minval(iu)
  ! the maximum/minimum value is used since a shifted (k+q)-mesh which is not
  ! commensurate can cause partially occupied states that are absent for the
  ! k-mesh
  iocc0=max(iocc0,iocc)
  iunocc0=min(iunocc0,iunocc)
  ! determine if system has a gap in energy
  ksgap=iocc0.lt.iunocc0
  if (iq.ne.0) then
     ! highest (partially) occupied state energy
     evlhpo=max(maxval(evalsv(iocc0,:)),maxval(evalsv0(iocc0,:)))
     ! lowest (partially) unoccupied state energy
     evllpu=min(minval(evalsv(iunocc0,:)),minval(evalsv0(iunocc0,:)))
  else
     ! highest (partially) occupied state energy
     evlhpo=maxval(evalsv(iocc0,:))
     ! lowest (partially) unoccupied state energy
     evllpu=minval(evalsv(iunocc0,:))
  end if
  ! check consistency with Fermi energy
  if (ksgap.and.((evlhpo.gt.efermi).or.(evllpu.lt.efermi))) then
     write(*,*)
     write(*,'("Error(findocclims): inconsistent Fermi energy for system&
          & with gap in energy:")')
     write(*,'(" Fermi energy            :",f12.6)') efermi
     write(*,'(" highest part. occ state :",f12.6)') evlhpo
     write(*,'(" lowest part. unocc state:",f12.6)') evllpu
     write(*,'(" recalculate Fermi energy or eigenvalues with reduced swidth")')
     write(*,*)
     call terminate
  end if
  ! *** assign nstval and nstcon ***
  nstval=max(iocc0,iocc)
  nstcon=nstsv-nstval
  if ((iocc0.ge.iunocc).or.(iocc.ge.iunocc0)) then
     write(*,'(a)') 'Info(findocclims): partially occupied states present'
  end if
  if (ksgap) then
     write(*,'(a)') 'Info(findocclims): system has KS-gap'
  else
     write(*,'(a)') 'Info(findocclims): no KS-gap found'
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
  if (.not.t) deallocate(evalsv0)
end subroutine findocclims
