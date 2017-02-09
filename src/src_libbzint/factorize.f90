!BOP
!
! !ROUTINE: factorize
!
! !INTERFACE:
      subroutine factorize(n,x,k,div)

! !DESCRIPTION: 
!
! This subroutine factorizes the real coordinates of a vector x, the output is an integer vector k, such that k(i)/div=x(i)
!
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: n

      real(8), intent(in) :: x(n)

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: div
      integer(4), intent(out) :: k(n)

! !DEFINDED PARAMETERS:

      integer(4), parameter :: nprim=168

! !LOCAL VARIABLES:

      integer(4) :: ik
      integer(4) :: ip
      integer(4), dimension(nprim) :: iprim

      integer(4) :: i,j
      integer(4) :: d(n)

      real(8) :: y,z
      logical :: found

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
!INTRINSIC ROUTINES:

! !EXTERNAL ROUTINES: 

! !REVISION HISTORY:
!       
! Created May 2006 by RGA

!EOP
!BOC

      do i=1,n
        found=.false.
        j=1
        do while (.not.found)
          j=j+1
          y=x(i)*dble(j)
          k(i)=int(y)
          z=dble(k(i))-y
          if(abs(z).lt.1.0d-10)then
            found=.true.
            d(i)=j
          endif
        enddo ! .not.found
      enddo ! i
      div=1
      do i=1,n
        div=div*d(i)
      enddo
      do i=1,n
        k(i)=k(i)*div/d(i)
      enddo

      do ip=1,nprim
        hec=div.ge.iprim(ip)
        do while (hec)
          hec = mod(div,iprim(ip)).eq.0
          ik=0
          do while((ik.lt.n).and.hec)
            ik=ik+1
            hec = hec.and.(mod(k(ik),iprim(ip)).eq.0)
          enddo
          if(hec) then
            div=div/iprim(ip)
            do ik=1,n
              k(ik)=k(ik)/iprim(ip)
            enddo
          endif
        enddo
      enddo
      if(div.eq.0) div=1
      return

      end subroutine factorize





