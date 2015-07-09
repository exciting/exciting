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
    integer(4) :: Recl
    real(8) :: eks, ehf, egw, sx, sc, vxc, deltax, deltae, z

! !DEFINED PARAMETERS:

    real(8), parameter :: tol=1.0d-3
 
! !EXTERNAL ROUTINES: 

    external writeqp

! !REVISION HISTORY:
!
! Created: 16.08.05 by RGA
! Readjusted: Mar 2013 by DIN
!
!EOP
!BOC      
    call cpu_time(tstart)
    if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
!
!   Allocate the array for the quasi-particle energies 
! 
    open(64,file='EVALQP.TXT',action='WRITE',form='FORMATTED')

    allocate(eqp(ibgw:nbgw,nkpt))

    do ikp = 1, nkpt
      write(64,1) ikp, vkl(:,ikp)
      write(64,2)
      do ie = ibgw, nbgw
!            
!       Calculate the new quasi-particle energy 
!
        sx=real(selfex(ie,ikp))
        sc=0.d0
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
               vxc = real(vxcnn(ie,ikp))+ ex_coef *real(selfex(ie,ikp))
           Else
               vxc = real(vxcnn(ie,ikp))
           End If
        Else         
           vxc = real(vxcnn(ie,ikp))
        End If 
        vxc=real(vxcnn(ie,ikp))
        z=0.d0        

        eks=evaldft(ie,ikp)

        deltax=sx-vxc
        ehf=eks+deltax

        egw=0.d0
        deltae=0.d0

        eqp(ie,ikp)=ehf

        write(64,3) ie, eks, ehf, egw, & 
        &           sx, sc, vxc, deltax, deltae, z        
        
      enddo ! ie
      write(64,*)
    enddo ! ikp
    close(64)

!----------------------------------------
!     Calculate QP Fermi energy to zero
!----------------------------------------
    nb=min(nbgw,int(chgval/2.d0)+30)
    call fermi(nkpt,nb-ibgw+1,eqp(ibgw:nb,:),ntet,tnodes,wtet,tvol, &
    &  nvelgw,.false.,eferqp,egap)

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
!     Repeat KS band structure analysis
!
      call bandanalysis('KS',ibgw,nbgw,evaldft(ibgw:nbgw,:),efermi)
!
!     QP band structure
!
      call bandanalysis('HF',ibgw,nbgw,eqp(ibgw:nbgw,:),eferqp)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      write(fgw,*)
      call write_cputime(fgw,tend-tstart, 'CALCEQPX')

    1 format('k-point #',i6,':',3f12.6)
    2 format(' state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk')    
    3 format(i4,'  ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)     

    return
end subroutine calceqpx
!EOC        
          
