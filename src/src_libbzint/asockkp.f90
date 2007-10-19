!BOP
!
! !ROUTINE: asockkp
!
! !INTERFACE:
      subroutine asockkp(wqk)
!
! !DESCRIPTION:
! Since the k-mesh and q-mesh are reduced by the symmetry. Without symmetry, k1
! and k2 are related by q1, k3 and k4 are related by q2. By using symmetry, k1 
! and k3 go into k, k2 and k4 go into k'. If both q1 and q2 go into q. Then we 
! will have some weight on wqk(q,k,k'). Since it represents at least two such kind
! of relations. wqk(q,k,k') is used to tell the wight of q on a specific pair of
! k and k'.This subroutine is only called by 'kqgen'. 
!       
! !USES:
      
      use kgen_internals     
     
      implicit none

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: wqk(:,:,:) ! the weight of each 
!                                             irreducible q-point 
!                                             for a given k and k' 
      
! !LOCAL VARIABLES:
      
      integer(4) :: i,j,m,ktmp(3),i1,j1,n1
      integer(4) :: idk1,idk2,idq
      integer(4), dimension(3) :: k1 !the k-point = 
      integer(4), dimension(3) :: k2 !the k'-point = k-q
      integer(4), dimension(3) :: q !the q-point 
      integer(4) :: nk               ! Total number of k-points
      integer(4), external :: idkp

!EOP
!BOC
      nk=div(1)*div(2)*div(3)
      wqk(:,:,:)=0.0d0
      do i=1,nirkp
        do j=1,nk
          do m=1,3
            k1(m)=ikp(m,i)
            k2(m)=kpt(m,j)
            ktmp(m)=k1(m)-k2(m)
            q(m)=mod(ktmp(m)+(1-isign(1,ktmp(m)))*div(m),div(m))
          enddo
          idq=idkp(q)
          idk1=idkp(k1)
          idk2=idkp(k2)
          i1=ikpid(redkp(idk1))
          j1=iqpid(redqp(idq))
          n1=ikpid(redkp(idk2))
          wqk(i1,n1,j1)=wqk(i1,n1,j1)+1
        enddo
      enddo
      
      end subroutine asockkp
!EOC      
          
        
