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
!      
!     Write quasi-particle energies to disk
!      
      call writeqpx

!     KS band structure
      call bandanalysis('KS',ibgw,nbgw,evaldft(ibgw:nbgw,:),efermi)

!     QP band structure
!     to calculate Fermi energy it is better to use 
!     only limited, low in energy, amount unoccupied states
      nb=min(nbgw,int(chgval/2.d0)+30)
      call fermi(nkpt,nb-ibgw+1,eqp(ibgw:nb,:),ntet,tnodes,wtet,tvol, &
      &  nvelgw,.false.,eferqp,egap)

      call bandanalysis('G0W0',ibgw,nbgw,eqp(ibgw:nbgw,:),eferqp)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'CALCEQP')

      return
      end subroutine calceqpx
!EOC        
          
