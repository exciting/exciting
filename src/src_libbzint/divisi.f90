!BOP
!
! !ROUTINE: divisi
!
! !INTERFACE:
      subroutine divisi(nkp,idiv,klist)

! !DESCRIPTION: 
!
! This subroutine factorizes the the integer coordinates of the k-points
! and their divisor (idiv) so that they are coprimes at output. That is,
! they have the minimum possible value      
!      
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: nkp ! Number of k-points
      
! !INPUT/OUTPUT PARAMETERS: 
     
      integer(4), intent(inout) :: idiv ! Minimum ommon divisor of the
!                                         integer k-points

      integer(4), intent(inout) :: klist(3,*) ! List of k-points

! !DEFINDED PARAMETERS:

      integer(4), parameter :: nprim=168

! !LOCAL VARIABLES:

      integer(4) :: ik
      integer(4) :: ir
      integer(4) :: ip
      integer(4), dimension(nprim) :: iprim

      logical hec


      data iprim /     2,  3,  5,  7, 11, 13, 17, 19, 23, 29, &
     &                31, 37, 41, 43, 47, 53, 59, 61, 67, 71, &
     &                73, 79, 83, 89, 97,101,103,107,109,113, &
     &               127,131,137,139,149,151,157,163,167,173, &
     &               179,181,191,193,197,199,211,223,227,229, &
     &               233,239,241,251,257,263,269,271,277,281, &
     &               283,293,307,311,313,317,331,337,347,349, &
     &               353,359,367,373,379,383,389,397,401,409, &
     &               419,421,431,433,439,443,449,457,461,463, &
     &               467,479,487,491,499,503,509,521,523,541, &
     &               547,557,563,569,571,577,587,593,599,601, &
     &               607,613,617,619,631,641,643,647,653,659, &
     &               661,673,677,683,691,701,709,719,727,733, &
     &               739,743,751,757,761,769,773,787,797,809, &
     &               811,821,823,827,829,839,853,857,859,863, &
     &               877,881,883,887,907,911,919,929,937,941, &
     &               947,953,967,971,977,983,991,997/

!EOP
!
!BOC
      do ip=1,nprim
        hec=idiv.ge.iprim(ip)
        do while (hec)
          hec = mod(idiv,iprim(ip)).eq.0
          ik=0
          do while((ik.lt.nkp).and.hec)
            ik=ik+1
            do ir=1,3
              hec = hec.and.(mod(klist(ir,ik),iprim(ip)).eq.0)
            enddo
          enddo
          if(hec) then
            idiv=idiv/iprim(ip)
            do ik=1,nkp
              do ir=1,3
                klist(ir,ik)=klist(ir,ik)/iprim(ip)
              enddo
            enddo
          endif
        enddo
      enddo
      if(idiv.eq.0) idiv=1
      end subroutine divisi
!EOC

