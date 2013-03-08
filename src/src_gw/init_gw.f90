!BOP
!
! !ROUTINE: init_gw
!
! !INTERFACE:
subroutine init_gw

! !DESCRIPTION:
!
! This is the main initialization subroutine of the gw program.
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
      integer(4) :: m, n
      logical    :: reduce
      real(8)    :: tstart, tend
      integer(4) :: stype_
 
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
      spinpol=associated(input%groundstate%spin)
      if (spinpol) then
        write(*,*)'GW EMERGENCY STOP!!!'
        stop 'Spin polarization is not implemented yet'
      end if

! Importantly, current implementation uses exclusively 
! "Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions", to be published in Comp. Phys. Commun. (2010)
      stype_=input%groundstate%stypenumber
      input%groundstate%stypenumber=-1

      if (.not.input%gw%skipgnd) then
!
! Generate eigenvectors for the complete (non-reduced) k-point set
!
        call system('cp INFO.OUT INFO.OUT-scf')
        reduce=input%groundstate%reducek
        input%groundstate%reducek=.false.          
        task=1
        input%groundstate%maxscl=1
        input%groundstate%nwrite=0
        call cpu_time(tstart)
        call gndstate
        call cpu_time(tend)
        if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
        call write_cputime(fgw,tend-tstart,'GNDSTATE')
        input%groundstate%reducek=reduce
      end if
      
! Generate the k- and q-point meshes
      call init_kqpts
      
!     if BSE is used just after GW, libbzint
!     may create the segmentation faults when try to use a general
!     k-point shift vkloff
!
      input%groundstate%stypenumber=stype_

! initialise the charge density and potentials from file
      Call readstate

! generate the core wavefunctions and densities
      Call gencore

! find the new linearisation energies
      Call linengy

! write out the linearisation energies
      Call writelinen

! generate the APW radial functions
      Call genapwfr

! generate the local-orbital radial functions
      Call genlofr
      
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
!     (still need to be checked!!!)
      call intipw

!     Generate the frequency mesh and integration weights.
      call init_freq

      return
end subroutine init_gw
!EOC
