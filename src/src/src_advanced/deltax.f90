
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!Created March 2008 (DIN)
!Modified Okt 2013 (DIN,UW) 
!Revision March 2014 (UW)

subroutine deltax
use modmain
use modmpi
implicit none
! local variables
integer is,ia,ias,ik
integer ir,irc,it,idm
integer :: msize
real(8) tau,resp,t1
! allocatable arrays
complex(8), allocatable :: vnlvv_full(:,:)

! Indexes of HOMO and LUMO orbitals
real(8), allocatable    :: delta(:,:)
integer    :: ihomo, ilumo,  iband
real(8)    :: k_Vxnl_k
complex(8) :: k_Vx_k
complex(8) zfinp
external   zfinp
! For matrix elements calculations
integer :: ilm
real(8), allocatable :: evalsvp(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
if (iscl.lt.1) return
!********************************************
! Calculation of non local matrix elements
!********************************************
allocate(vnlvv_full(nstsv,nkpt))
#ifdef MPI
         Do ik = firstk (rank), lastk (rank)
#else
         Do ik = 1, nkpt
#endif       
             call oepvnlk_deltax(ik,vnlvv_full(1,ik))     
         End Do
#ifdef MPI
         msize = nstsv
       call mpi_allgatherv_ifc(nkpt,msize,zbuf=vnlvv_full)
#endif

!********************************************
! Calculation of the potential discontinuity
!         after last iteration
!********************************************
  ! For matrix elements calculations
  allocate(evalsvp(nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(evecfv(nmatmax,nstfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngrtot,nspinor,nstsv))
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
  allocate(zrhoir(ngrtot))
  allocate(delta(nstsv,nkpt))
  !write(*,*) "EFERMI = ", efermi
  !write(*,*) "nstsv = ", nstsv

#ifdef MPI
Do ik = firstk (rank), lastk (rank)
#else
Do ik = 1, nkpt
#endif       

     call getevalsv(vkl(1,ik),evalsvp)
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call getevecsv(vkl(1,ik),evecsv)
     ! find the matching coefficients
     call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
     ! calculate the wavefunctions for all states for the input k-point
     call genwfsv(.false.,ngk(ik,1),igkig(1,ik,1),evalsvp,apwalm,evecfv,evecsv,& 
      wfmt,wfir)

!    For each band calculate <k|Vx_nl - Vx|k>
     do iband=1,nstsv

        k_Vxnl_k=0d0
        k_Vx_k=0d0

    
        ! 1 step:
        ! <k|Vx_nl|k>
        ! The calculations are done in oepvnkl.f90
        ! Output: k_Vxnl_k(nkpt)
        k_Vxnl_k=vnlvv_full(iband,ik)

        ! 2 step:
        ! Calculation of  <k|Vx|k>
        ! get the eigenvalues/vectors from file for input k-point

        !  ------ <k|Vx|k> ------
        do is=1,nspecies
         do ia=1,natoms(is)
            ias=idxas(ia,is)
            do ilm=1,lmmaxvr
            do irc=1,nrcmtmax
               zrhomt(ilm,irc,ias)=zvxmt(ilm,irc,ias)*wfmt(ilm,irc,ias,1,iband)
            end do
            end do
         end do
        end do
        do ir=1,ngrtot
           zrhoir(ir)=zvxir(ir)*wfir(ir,1,iband)
        end do
        k_Vx_k=zfinp(.false.,wfmt(1,1,1,1,iband),zrhomt,&
                             wfir(1,1,iband),zrhoir)

        delta(iband,ik) = (k_Vxnl_k-real(k_Vx_k))

      end do ! iband

!  end do
   End Do  ! kpoints
   
#ifdef MPI
       call mpi_allgatherv_ifc(nkpt,nstsv,rbuf=delta)
#endif
! Write delta into file
If (rank .Eq. 0) then 
   open(500,file='DELTAX'//trim(filext),action='WRITE',form='FORMATTED')
   write(500,*) nkpt, nstsv
   Do ik = 1, nkpt
      write(500,'(3f12.8,i5)') vkl(1,ik), vkl(2,ik), vkl(3,ik)
      do iband=1,nstsv
        write(500,'(6e20.6)') delta(iband,ik)
      end do
   end Do
  close(500)
End if
  deallocate(evalsvp,evecfv,evecsv)
  deallocate(apwalm)
  deallocate(wfmt,wfir)
  deallocate(zrhomt,zrhoir)
  deallocate(delta)
!end if

!********************************************
deallocate(vnlvv_full)
return
end subroutine

