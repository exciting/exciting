!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: mpiresumeevec
! !INTERFACE:
!
!
Subroutine mpiresumeevecfiles
      Use modmain
#ifdef MPI
      Use modmpi
! !DESCRIPTION:
!  Routine reads the lokal EIGVECk1-k2.OUT files that were created to
!  avoid file system inconsistencies. And writes them into the standard
!  EIGVEC.OUT File
!
! !REVISION HISTORY:
!   Created October SEPT 2006 (MULEOBEN)
!   by Cristian Meisenbichler
!EOP
      Implicit None
  ! arguments
      Integer :: ik, proc, nmatmax_, nstfv_, nspnfv_, nstsv_, recl
      Integer :: recvstatus (MPI_STATUS_SIZE),neighbors(procs)
      Character (256) :: filetag
      Complex (8) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8) :: evecsv (nstsv, nstsv)
      Real (8) :: evalfv (nstfv, nspnfv), vkl_ (3), evalsvp (nstsv), &
     & occsvp (nstsv)
      Character (256), External :: outfilenamestring

      If (splittfile .And. (procs .Gt. 1) ) Then
! start a receive in order to pass around a token from rank 0 to max

         filetag = 'EVECFV'
         Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, &
        & evecfv
         call resumefile(filetag,Recl)

         filetag = 'EVECSV'
         Inquire (IoLength=Recl) vkl_, nstsv_, evecsv
         call resumefile(filetag,Recl)

         filetag = 'EVALFV'
         Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_, evalfv
         call resumefile(filetag,Recl)

!
         filetag = 'EVALSV'
         Inquire (IoLength=Recl) vkl_, nstsv_, evalsvp
         call resumefile(filetag,Recl)

!
         filetag = 'OCCSV'
         Inquire (IoLength=Recl) vkl_, nstsv_, occsvp
         call resumefile(filetag,Recl)

!
         Call SYSTEM ("sync")
         !If (rank==0) Then
         !   Write (60,*) "resumed split files"
         !   Call flushifc (60)
         !End If

         !if I am not the last process pass on the token
      End If
      call barrier
      Call SYSTEM ("sync")
      splittfile = .False.

#endif
!
      Return
End Subroutine mpiresumeevecfiles

subroutine resumefile(filetagarg,Recl)
use modmpi
use mod_misc
use mod_kpoint
use modinput

	 Implicit None
	 !arguments
	 integer,parameter::filetaglenth=256
	character(filetaglenth)::filetagarg
	integer::Recl,reclloc,recordunit_inbytes
	character,allocatable::buffer(:)
	!local
	integer::ik
#ifdef MPI
	Character (256), External :: outfilenamestring

	!compute record unit in bytes abusing the filetag variable
	Inquire (IoLength=Reclloc) filetagarg
	recordunit_inbytes=filetaglenth/Reclloc

allocate(buffer(Recl*recordunit_inbytes))

    if(rank.eq.0 .or. (.not. input%sharedfs .and. firstinnode)) then
     Open (71, File=trim(filetagarg)//trim(filext), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
	endif
    if(rank.lt.nkpt)Open (77, File=outfilenamestring(filetagarg, firstk(rank)), &
        & Action='READ', Form='UNFORMATTED', Access='DIRECT', &
        & Recl=Recl)
    Do ik =1,nkpt
    	if(procofk(ik).eq.rank .and. (rank.lt.nkpt)) then
            Read (77, Rec=ik-firstk(procofk(ik))+1) buffer
        endif
        Call MPI_bcast (buffer,Recl*recordunit_inbytes , MPI_CHARACTER, procofk(ik), MPI_COMM_WORLD, ierr)

        if(rank.eq.0 .or. (.not. input%sharedfs .and. firstinnode)) then
            Write (71, Rec=ik) buffer
		endif
    End Do
    Close (77, Status='DELETE')
    Close (71)
    deallocate(buffer)
#endif
end subroutine
