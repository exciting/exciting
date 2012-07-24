!BOP
!
! !ROUTINE: initgw
!
! !INTERFACE
      subroutine initgw

! !DESCRIPTION:
!
! This is the main initialization subroutine of the gw program,
! 
! !USES:
      use modmain
      use modgw

! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: i
      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: ic
      integer(4) :: il
      integer(4) :: is      
      integer(4) :: ist      
      integer(4) :: m
 
! !EXTERNAL ROUTINES: 

      external init0
      external init1
      external readstate
      external genvxcig
      external init_freq
      external readingw
      external initkqpts

! !REVISION HISTORY:
!
! Created 16. May. 2006 by RGA
! Revisited, DIN: 26.04.2011
!
!EOP
!BOC
      if (.not.input%groundstate%tetra) then
        write(*,*)'GW EMERGENCY STOP!!!'
        write(*,*)'k-point meshes should be generated with tetra =.true.'
        stop 'ERROR in initgw'
      endif
      spinpol=associated(input%groundstate%spin)
      if (spinpol) then
        write(*,*)'GW EMERGENCY STOP!!!'
        stop 'Spin polarization is not implemented yet'
      end if

!     Minimal muffin-tin radius
      rmtmin=rmt(1)
      if(nspecies.gt.1)then
        do is=2,nspecies
          if (rmt(is).lt.rmtmin) rmtmin=rmt(is)
        enddo
      endif

!     Determine G-vector cutoff parameters
      gkmax=input%groundstate%rgkmax/rmtmin
      gqmax=input%gw%MixBasis%gmb*gkmax
      gmaxbarc=min(pwm*gqmax,input%groundstate%gmaxvr)
      
      if(gmaxbarc.gt.input%groundstate%gmaxvr)then
        write(*,*)'WARNING(initgw)! One should increase the value of gmaxvr:'
        write(*,*) 'gkmax=',gkmax,'    gqmax=', gqmax
        write(*,*) 'gmaxvr',input%groundstate%gmaxvr,'    gmaxbarc=', gmaxbarc
      end if
      
!     initialise universal variables
      call init0
      call init1

!     initialize or read the charge density and potentials from file
      call readstate

!     Generate the frequency mesh and integration weights.
      call init_freq
      
!     Generate the k- and q-point meshes
      call initkqpts

!     Tranform xc potential to reciprocal space
      call genvxcig

!     determine the number of core states for each species (auxiliary arrays)
      allocate(ncore(nspecies))
      ncmax=0
      nclm=0
      ncg=0
      lcoremax=0
      do is=1,nspecies
        ncore(is)=0
        ic = 0
        do ist=1,spnst(is)
          if (spcore(ist,is)) then
            ncore(is)=ncore(is)+1
            il=spl(ist,is)
            do m=-spk(ist,is),spk(ist,is)-1
              ic=ic+1
            end do
          end if
        end do
        ncmax=max(ncmax,ncore(is))
        nclm=max(nclm,ic)
        lcoremax=max(lcoremax,il)
        ncg=ncg+ic*natoms(is)
      end do

!     setting a unique index for all the core states of all atoms
      call setcorind

!     reciprocal cell volume
      vi=1.0d0/omega

!     shortcut for basis vectors 
      avec(:,1)=input%structure%crystal%basevect(:,1)
      avec(:,2)=input%structure%crystal%basevect(:,2)
      avec(:,3)=input%structure%crystal%basevect(:,3)
      
!     reciprocal lattice basis lengths
      do i=1,3
         alat(i)=dsqrt(avec(1,i)*avec(1,i)+ &
                       avec(2,i)*avec(2,i)+ &
                       avec(3,i)*avec(3,i))
         pia(i)=2.0d0*pi/alat(i)
      end do

!     additional arrays used for convenience
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
!         shortcut for atomic positions
          atposl(:,ia,is)=input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
        end do
!       calculate the muffin-tin volume
        vmt(is)=4.0d0*pi*rmt(is)*rmt(is)*rmt(is)/(3.0d0*omega)
      end do

!     Calculate the overlap between two PW 
!     In exciting there is the same quantity: conjg(cfunig(:)) = ipwint(:)
!     But beed to be checked!!!
      call intipw

      return
      end subroutine initgw
!EOC
