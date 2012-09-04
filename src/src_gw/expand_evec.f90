!BOP
!
! !ROUTINE: expand_evec
!
! !INTERFACE:
      subroutine expand_evec(ik,trans)

! !DESCRIPTION:

!Calculates the coeficients of the expansion of the eigenvectors of the DFT 
!calculation within the MT-spheres in the spherical basis as follows:
!
!\begin{equation}
!\Psi_{n\vec{k}}(\vec{r})=\sum\limits_a{\sum\limits_{l=0}^{lmax}{%
!\sum\limits_{m=-l}^{l}{\left[\mathcal{A}_{lm}^{na}(\vec{k})u_{al}(r)+%
!\mathcal{B}_{lm}^{na}(\vec{k})\dot{u}_{al}(r)+\mathcal{C}_{lm}^{na}(\vec{k})u_{al}(r,E_2)%
!\right]Y_{lm}(\hat{r})}}}
!\end{equation}
!where $a$ indicates the atoms. The corresponding expresions for the
!coefficients are:
!
!\begin{subequations}\label{almevec}
!\begin{align}
!\mathcal{A}_{lm}^{na}(\vec{k})=&\sum\limits_{\vec{G}}{Z^n_{\vec{k}+\vec{G}}%
!A^a_{lm}(\vec{k}+\vec{G})}\\
!\mathcal{B}_{lm}^{na}(\vec{k})=&\sum\limits_{\vec{G}}{Z^n_{\vec{k}+\vec{G}}%
!B^a_{lm}(\vec{k}+\vec{G})}\\
!\mathcal{C}_{lm}^{na}(\vec{k})=&Z^n_{\vec{k}+\vec{G}_{LO}}%
!C^a_{lm}(\vec{k}+\vec{G}_{LO})
!\end{align}
!\end{subequations}
!
! The respective coefficients are stored as:
!
!\begin{subequations}
!\begin{align}
!\mathcal{A}_{lm}^{na}(\vec{k})&\rightarrow\texttt{alfa(a,n,i)}\\
!\mathcal{B}_{lm}^{na}(\vec{k})&\rightarrow\texttt{beta(a,n,i)}\\
!\mathcal{C}_{lm}^{na}(\vec{k})&\rightarrow\texttt{gama(a,n,i)}
!\end{align}
!\end{subequations}
! 
! with $i=l^2+l+m+1$
!
! !USES:
      
      use modinput
      use modmain
      use modgw      

! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: ik    ! index of the k-point      
      character(1), intent(in) :: trans
      
! !LOCAL VARIABLES:

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: io
      integer(4) :: is
      integer(4) :: ist
      integer(4) :: l1    ! orbital angular momentum quantum number
      integer(4) :: l1m1  ! l1^2+l1+m1+1
      integer(4) :: m1    ! z component of the orbital ang. mom
     
      complex(8), allocatable :: apwalm(:,:,:,:,:)
 
! !EXTERNAL ROUTINES: 

      complex(8), external :: zdotc
      complex(8), external :: zdotu
      external match

! !REVISION HISTORY:
! 
! Created Dic. 2003 by RGA
! Last Modification: July. 20th. 2004 by RGA
! Revision: 5.05.2011 by DIN
!
!EOP  
!BOC  
!     allocate the arrays to store the A_lm, B_lm and C_lm coefficients
      allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))

!     find the matching coefficients
      call match(ngknr(1,ik),gkcnr(:,1,ik),tpgkcnr(:,:,1,ik), &
     &  sfacgknr(:,:,1,ik),apwalm(:,:,:,:,1))

      select case (trans)
      case ('t','T')
        do is=1,nspecies
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            do l1=0,input%groundstate%lmaxapw
              do m1=-l1,l1
                l1m1=idxlm(l1,m1)
                do io=1,apword(l1,is)
                  do ist=1,nstfv
                    eveckalm(ist,io,l1m1,ias,1)=zdotu(ngknr(1,ik), &
     &                  eveck(1,ist,1),1,apwalm(1,io,l1m1,ias,1),1)
                  enddo ! ist  
                enddo ! io
              enddo ! m1     
            enddo ! l1       
          enddo  !ia
        enddo !is
      case ('c','C')
        do is=1,nspecies
          do ia = 1, natoms(is)
            ias=idxas(ia,is)
            do l1=0,input%groundstate%lmaxapw
              do m1=-l1,l1
                l1m1=idxlm(l1,m1)
                do io=1,apword(l1,is)
                  do ist=1,nstfv
                    eveckpalm(ist,io,l1m1,ias,1)=zdotc(ngknr(1,ik), &
     &                  apwalm(1,io,l1m1,ias,1),1,eveckp(1,ist,1),1)
                  enddo ! ist  
                enddo ! io
              enddo ! m1     
            enddo ! l1       
          enddo  !ia
        enddo !is
      case default
        write(*,*)'ERROR in expand_evec'
        write(*,*)'Allowed values of trans are t, T, c or C'
        write(*,'(a19,a1)')'Received trans = ',trans 
        stop 
      end select
      deallocate(apwalm)

      return
      end subroutine expand_evec
!EOC      
