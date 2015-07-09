!BOP
! !ROUTINE: inimb
!
! !INTERFACE:
      subroutine init_mb

! !DESCRIPTION:
!
! This subroutine initializes the mixed basis. Calculates the product
! functions and their overlap matrix, and orthonormalizes them. Calculates
! the matrices rlam, rrint , bradket and uxcu needed afterwards.
! 
! !USES:

      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: irm
      integer(4) :: is
      integer(4) :: lms
      integer(4) :: ist, l, m, lm, ir, im
      real(8)    :: norm
      logical    :: core_ortho

! !EXTERNAL ROUTINES: 
      
      external gencore
      external linengy
      external writelinen
      external genapwfr
      external genlofr
      external diagsgi
      external setumix
      external calcrlamint
      external calcbradket
      external calctildeg

! !INTRINSIC ROUTINES: 

      intrinsic log
      intrinsic exp
      intrinsic atan

! !REVISION HISTORY:
!
! Last Modified 10.11.2005 by RGA
! Revisited 28.04.2011 by DIN
!      
!EOP
!BOC
      if (debug) open(701,file='MIXEDBASIS.OUT',action='WRITE',status='UNKNOWN',form='FORMATTED')

!------------------------------------------------------------------------------

      if (iopcore<3) then
        ! According to the definition of core wafefunction in FHIgap code [Eq.(1.1.3)],
        ! one has to include the following prefactor into radial part.
        ! In addition I change the EXCITING definition, where ucore = r*rwfcr
        if (allocated(ucore)) deallocate(ucore)
        allocate(ucore(spnrmax,2,spnstmax,natmtot))
        do is=1,nspecies
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            do ist=1,ncore(is)
              l=spl(ist,is)
              norm=sqrt(0.5d0*spocc(ist,is)/dble(2*l+1))
              do ir=1,nrmt(is)
                ucore(ir,1,ist,ias)=norm*rwfcr(ir,1,ist,ias)/spr(ir,is)
              end do
            enddo ! ist
          end do
        end do
        !core_ortho=.false.    
        !if(core_ortho) call orthog_corewf
      end if
      
!------------------------------------------------------------------------------

!     Generate all possible radial function products
      call setuprod

!     calculate the radial part of the mixed basis functions
      call setumix
      
!     Calculate the matrix elements <r^L>_{aNL} and
!     <r_<^L/r_>^{L+1}> _{aNL,N'L}:
      call calcrlamint

!     Calculate the matrix elements <NL,lambda|lambda'>:
      call calcbradket
      
!     Calculate the total number of mixed wave functions (including M)    
!     = size of the local part of the matrices
      locmatsiz=0
      lmixmax=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          lms=0
          do irm=1,nmix(ias)
            lms=lms+(2*bigl(ias,irm)+1)
          enddo  
          if (lms.gt.lmixmax) lmixmax=lms
          locmatsiz=locmatsiz+lms
        enddo ! ia
      enddo ! ias

!     The maximum size of MB basis
      matsizmax=locmatsiz+ngqmax   
      
      call setlocmixind
      
      if (maxbigl>input%groundstate%lmaxapw) then
        if (allocated(idxlm)) deallocate (idxlm)
        allocate(idxlm(0:maxbigl,-maxbigl:maxbigl))
        lm = 0
        do l = 0, maxbigl
          do m = -l, l
            lm = lm + 1
            idxlm(l,m) = lm
          end do
        end do
      end if
      
      !-------------------------------------------------------------------
      ! mapping: MB function index -> (aNLM)
      !-------------------------------------------------------------------
      if (allocated(mbindex)) deallocate(mbindex)
      allocate(mbindex(locmatsiz,5))
      mbindex(:,:) = 0
      im = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do irm = 1, nmix(ias)
            l = bigl(ias,irm)
            do m = -l, l
              im = im+1
              mbindex(im,1) = is
              mbindex(im,2) = ia
              mbindex(im,3) = irm
              mbindex(im,4) = l
              mbindex(im,5) = m
            end do ! m
          end do ! irm
        end do ! ia
      end do ! is

!     Calculate the coefficients tildeg needed for the structure constants
      call calctildeg(2*(input%gw%MixBasis%lmaxmb+1))
      
      if(debug)close(701)
      
      return
      end subroutine
!EOC
