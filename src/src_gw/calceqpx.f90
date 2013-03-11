!BOP
!
! !ROUTINE: calceqpx
!
! !INTERFACE:
      subroutine calceqpx

! !DESCRIPTION:
! 
! Given the matrix elements $\langle
! \Psi_{n\vec{k}}|\Sigma(\vec{k},\omega)|\Psi_{n\vec{k}}\rangle$, $\langle
! \Psi_{n\vec{k}}|V^{xc}|\Psi_{n\vec{k}}\rangle$ and
! $\varepsilon^{DFT}_{n\vec{k}}$. this subroutine calculates the
! quasi-particle energies $\varepsilon^{qp}_{n\vec{k}}$
!
! !USES:
!
      use modmain
      use modgw     
       
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: ie   !(Counter) Runs over bands
      integer(4) :: ikp  !(Counter) Runs over k-points
      integer(4) :: nb
      real(8) :: delta   ! absolute value of the difference between
!                          quasiparticle energies of succesive iterations 
      real(8) :: eferqp, egap
      real(8) :: tstart, tend
      integer(8) :: Recl

! !DEFINED PARAMETERS:

      real(8), parameter :: tol=1.0d-3
 
! !EXTERNAL ROUTINES: 

      external writeqp

! !REVISION HISTORY:
!
! Created: 16.08.05 by RGA
!
!EOP
!BOC      
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
!
!     Allocate the array for the quasi-particle energies 
! 
      allocate(eqp(ibgw:nbgw,nkpt))
!        
!     Loop over kpoints
! 
      do ikp = 1, nkpt
!         
!       Loop over bands
! 
        do ie = ibgw, nbgw
!            
!         Calculate the new quasi-particle energy 
!
          delta=real(selfex(ie,ikp))-real(vxcnn(ie,ikp))
          eqp(ie,ikp)=evaldft(ie,ikp)+delta
        enddo ! ie
      enddo ! ikp

!     QP band structure
!     to calculate Fermi energy it is better to use 
!     only limited, low in energy, amount unoccupied states
      nb=min(nbgw,int(chgval/2.d0)+30)
      call fermi(nkpt,nb-ibgw+1,eqp(ibgw:nb,:),ntet,tnodes,wtet,tvol, &
      &  nvelgw,.false.,eferqp,egap)

!----------------------------------------
!     Set QP fermi energy to zero
!----------------------------------------
      eqp(:,:)=eqp(:,:)-eferqp
      eferqp=0.d0

!----------------------------------------
!     Save QP energies into binary file
!----------------------------------------
      Inquire (IoLength=Recl) nkpt, ibgw, nbgw, &
     &  vkl(:,1), eqp(ibgw:nbgw,1), evaldft(ibgw:nbgw,1)
      Open (70, File='EVALQP.OUT', Action='WRITE', Form='UNFORMATTED', &
     &   Access='DIRECT', status='REPLACE', Recl=Recl)
      do ikp = 1, nkpt
        write(70, Rec=ikp) nkpt, ibgw, nbgw, &
     &    vkl(:,ikp), eqp(ibgw:nbgw,ikp), evaldft(ibgw:nbgw,ikp)
      end do ! ikp
      Close(70)
!      
!     Write quasi-particle energies into QPENE-eV.OUT
!      
      call writeqpx
!      
!     Repeat KS band structure analysis
!
      call bandanalysis('KS',ibgw,nbgw,evaldft(ibgw:nbgw,:),efermi)
!
!     QP band structure
!
      call bandanalysis('G0W0',ibgw,nbgw,eqp(ibgw:nbgw,:),eferqp)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      write(fgw,*)
      call write_cputime(fgw,tend-tstart, 'CALCEQP')

      return
      end subroutine calceqpx
!EOC        
          
