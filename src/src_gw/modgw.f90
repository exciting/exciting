!BOP
! !MODULE: modgw
! !DESCRIPTION:
!   Contains all the global variables required by the GW EXCITING code.
!
! !REVISION HISTORY:
!
!   Created MAY 2006 (RGA)
!   Revisited, 27.04.2011 by DIN
!
#include "maxdefinitions.inc"
module modgw

#ifdef TETRA
      use modtetra
#endif

! !PUBLIC VARIABLES:

!----------------------------!
!     general                !
!----------------------------!
! Task identification number
      integer(4) :: testid
! general gw output file
      integer(4) :: fgw
! shortcut for atomic position array
      real(8) :: atposl(3,_MAXATOMS_,_MAXSPECIES_)
! shortcut for basis vectors
      real(8) :: avec(3,3)
! lengths of the basis vectors
      real(8) :: alat(3)      
! 2*pi/a, 2*pi/b, 2*pi/c
      real(8) :: pia(3)
! 1/omega - reciprocal unitcell volume
      real(8) :: vi 
! for more detailed (debug) information
      logical :: debug
! muffin-tin volume, relative to cell volume
      real(8) :: vmt(_MAXSPECIES_)
! spin-polarized flag
      logical :: spinpol
! spin factor to treat the spin degeneracy
      real(8) :: sfact=2.0d0
! Ha --> eV
      real(8), parameter :: hev=27.2113961

!----------------------------!
!     frequency grid         !
!----------------------------!
! Number of frequencies
      integer(4) :: nomeg 
! Flag for the type of frequency mesh
      integer(4) :: wflag
! cutoff on frequency axis
      real(8) :: freqmax 
! frequency grid
      real(8), allocatable :: freqs(:)
! integration weights for frequency grid
      real(8), allocatable :: womeg(:) 

!--------------------------------------------!
!     Brillouin zone integration             !
!--------------------------------------------!
! flag for the frequency dependence of the  q dependent integration.
!              =1 No frequency dependence included
!              =2 Real Frequencies
!              =3 Imaginary frequencies
      integer(4) :: fflg
! flag for the q dependent integration.
!              = 0  direct sum over k-points
!              = 1 linearized tetrahedron method
      integer(4) :: convflg
! Weights for integration over unoccupied states
      real(8), allocatable :: unw(:,:,:,:,:)   
! Weights for convolutions when k is the integration variable
      real(8), allocatable :: kcw(:,:,:,:) 
! Weights for surface integration over unoccupied states
      real(8), allocatable :: unwsurf(:,:,:,:,:) 
! Weights for surface convolutions when k is the integration variable
      real(8), allocatable :: kcwsurf(:,:,:,:) 

      real(8), allocatable :: kiw(:,:)
      real(8), allocatable :: kwfer(:,:)
      real(8), allocatable :: ciw(:,:)
      
      logical :: metallic
      logical :: Gamma

!-----------------------------------!
!     tetrahedron method variables  !
!-----------------------------------!
! the commented variables are defined in modtetra.f90
#ifndef TETRA
! integer k-point offset
      integer(4) :: ikloff(3)
! k-point offset divisor
      integer(4) :: dkloff
! k-points common divisor
      integer(4) :: dvk
! Number of tetrahedra
      integer(4) :: ntet
! index of the k-points corresponding to the nodes of each tetrahedron
      integer(4), allocatable :: tnodes(:,:)
! weight of each tetrahedron.
      integer(4), allocatable :: wtet(:)
! volume of the tetrahedra relative to the BZ volume
      real(8)    :: tvol      
      integer(4) :: mnd
! k-dependent weight of each q-point
      integer(4), allocatable :: kqid(:,:)
! q-points common divisor
      integer(4) dvq
! auxiliary arrays
#endif
      
      integer(4), allocatable :: indkp(:)
      integer(4), allocatable :: idikp(:)     ! kpoint index of the irred. kpoint

      integer(4) :: ntetnr                    ! Total number of tetrahedra
      integer(4), allocatable :: wtetnr(:)    ! weight of each tetrahedron  for integration
      integer(4), allocatable :: tnodesnr(:,:)! index of the k-points corresponding to the nodes of each tetrahedra for integration

      integer(4), allocatable :: linkq(:,:)
      integer(4), allocatable :: iwkp(:)

!------------------------------!
!     Non-reduced G+k arrays   !
!------------------------------!
! number of G+k-vectors for augmented plane waves
      Integer, Allocatable :: ngknr(:,:)
! index from G+k-vectors to G-vectors
      Integer, Allocatable :: igkignr(:,:,:)
! G+k-vectors in lattice coordinates
      Real (8), Allocatable :: vgklnr(:,:,:,:)
! G+k-vectors in Cartesian coordinates
      Real (8), Allocatable :: vgkcnr(:,:,:,:)
! length of G+k-vectors
      Real (8), Allocatable :: gkcnr(:,:,:)
! (theta, phi) coordinates of G+k-vectors
      Real (8), Allocatable :: tpgkcnr(:,:,:,:)
! structure factor for the G+k-vectors
      Complex (8), Allocatable :: sfacgknr(:,:,:,:)

!--------------------------------!
!     Small group of q-vectors   !
!--------------------------------! 
!     non-reduced number of q-points
      integer :: nqptnr
!     number of the symmetry operations in the small group of q
      integer, allocatable :: nsymq(:) 
!     q-dependent k-point weight
      real(8), allocatable :: wkpq(:,:)

!     number of k-points in IBZ(q)
      integer, allocatable :: nkptq(:)
!     index of the symmetry operation which rotates the k-point into equivalent one
      integer, allocatable :: iksymq(:,:)
!     map the k-point index to the corresponding irreducible one
      integer, allocatable :: indkpq(:,:)
!     map the irreducible k-point index to the corresponding from the non-reduced set
      integer, allocatable :: idikpq(:,:)
!     rotation matrix for ylm's      
      complex(8), allocatable :: djmm(:,:)
      
      integer, allocatable :: nsymkstar(:,:), isymkstar(:,:,:)
      integer, allocatable :: n12dgn(:,:,:)
      
!----------------------------!
!     LAPW eigenfunctions    !
!----------------------------!
! LAPW eigenenergies
      real(8), allocatable :: evaldft(:,:)
! Eigenvectors at k
      complex(8), allocatable :: eveck(:,:,:)
! Eigenvectors at k'=k-q
      complex(8), allocatable :: eveckp(:,:,:)
! Spherical harmonic expansion coefficients of the eienvector at k
      complex(8), allocatable :: eveckalm(:,:,:,:,:)      
! Spherical harmonic expansion coefficients of the eienvector at k'=k-q
      complex(8), allocatable :: eveckpalm(:,:,:,:,:)      
! Highest partially occupied band      
      integer(4) :: nomax
! Lowest partially unoccupied band      
      integer(4) :: numin           
! Lower band index for GW output
      integer(4) :: ibgw
! Upper band index for GW output
      integer(4) :: nbgw
! Number of bands for gw output      
      integer(4) :: nbandsgw
! Number of electrons used in GW
      real(8)    :: nvelgw   

!----------------------------!
!     Core states            !
!----------------------------!
      integer :: iopcore   ! option for core treatment 
                           !  iopcore  = 0  --- core states are included in all calculations 
                           !           = 1  --- core states are included in exchange  but not in correlation 
                           !           = 2  --- core states are excluded in all calculations, but kept in 
                           !                    the construction of mixed basis 
                           !           = 3  -- core states are excluded completely
! maximum number of core states per atom
      integer(4) :: ncmax
! Max. num of core states including lm
      integer(4) :: nclm 
! Max. L of core states
      integer(4) :: lcoremax      
! Total num of core states over all atoms including lm
      integer(4) :: ncg  
! number of core states of each species
      integer(4), allocatable :: ncore(:)  
! indexes of the core states
!     corind(:,1)= species
!     corind(:,2)= atom of species
!     corind(:,3)=core state
!     corind(:,4)=l
!     corind(:,5)=m
      integer(4), allocatable :: corind(:,:)
      real(8), allocatable :: ucore(:,:,:,:)

!----------------------------!
!     Local orbitals         !
!----------------------------!
      real(8), allocatable :: lorwf(:,:,:,:,:)
        
!----------------------------!
!     XC potential           !
!----------------------------!
! G-space interstitial exchange-correlation potential
      complex(8), allocatable :: vxcig(:)
! APW-APW exchange-correlation  integrals
      real(8), allocatable :: vxcraa(:,:,:,:,:,:)
! local-orbital-APW exchange-correlation  integrals
      real(8), allocatable :: vxcrloa(:,:,:,:,:)
! local-orbital-local-orbital exchange-correlation  integrals
      real(8), allocatable :: vxcrlolo(:,:,:,:)
! diagonal matrix elements of the exchange-correlation potential
      complex(8), allocatable :: vxcnn(:,:)

!----------------------------!
!     mixed basis (general)  !
!----------------------------!
! Size of the mixed basis
      integer(4) :: matsiz, matsizmax
! Matrix elements M^i_nm and \tilde{M}^i_nm
      complex(8), allocatable :: minmmat(:,:,:)
! Matrix elements M^i_cm and \tilde{M}^i_cm
      complex(8), allocatable :: micmmat(:,:,:)
! Matrix elements M^i_nc and \tilde{M}^i_nc
      complex(8), allocatable :: mincmat(:,:,:)
! Matrix elements M^i_cc and \tilde{M}^i_cc
      complex(8), allocatable :: miccmat(:,:,:)
! Matrix elements between mixed functions and constant function
      complex(8), allocatable :: wi0(:)
! Matrix elements between mixed functions and planewaves
      complex(8), allocatable :: mpwmix(:,:)      
! Overlap of two PW
      complex(8), allocatable :: ipwint(:)

!----------------------------!
!     mixed basis (MT Spheres)!
!----------------------------!
! Maximum L of the left product functions
      integer(4) :: lleftmax
! Radial product functions
      real(8), allocatable :: uprod(:,:,:)
! Overlap matrix of the product functions
      real(8), allocatable :: umat(:,:,:)
! Maximum number of product functions
      integer(4) :: maxnup
! Number of product functions per atom
      integer(4), allocatable :: nup(:)
! l,l' pairs of the product functions
      integer(4), allocatable :: eles(:,:,:)
! Radial mixed functions
      real(8), pointer :: umix(:,:,:)
! L quantum number of the mixed functions
      integer(4), pointer :: bigl(:,:)
! Maximum L of the mixed functions per atom
      integer(4), allocatable :: mbl(:)
! Maximum L of the mixed functions
      integer(4) :: maxbigl
! Number of mixed radial functions per atom
      integer(4), allocatable :: nmix(:)
! Maximum number of radial functions per atom
      integer(4) :: maxnmix
! <umix(l)|r^(l+2)> integrals     
      real(8), allocatable :: rtl(:,:)
! <umix(l1)|r^(l1)/r^(l2+1)|umix(l2)> integrals
      real(8), allocatable :: rrint(:,:)
! <umix(L)|ucore(l1)u(l2)> integrals
      real(8), allocatable :: bradketc(:,:,:,:,:,:)
! <umix(L)|u(l1)u(l2)> integrals
      real(8), allocatable :: bradketa(:,:,:,:,:,:,:)
! <umix(L)|ulo(l1)u(l2)> integrals
      real(8), allocatable :: bradketlo(:,:,:,:,:,:)
! the gaunt coefficients
      real(8), allocatable :: cgcoef(:) 
! Size of the local part of the mixed basis
      integer(4) :: locmatsiz
! maximun number of mixed functions per atom
      integer(4) :: lmixmax 
! indexes of the local mixed basis functions
      integer(4), allocatable :: locmixind(:,:)
      integer(4), allocatable :: mbindex(:,:)
! The matrix elements jlam      
      real(8), allocatable  :: jlam(:,:) 

!---------------------------------!
!     mixed basis (Interstitial)  !
!---------------------------------!
! G-vector cut-off for the mixed basis
      real(8) :: gqmax
!  number of G+q-vectors for the mixed basis
      integer(4), allocatable :: ngq(:)
! maximum number of G+q-vectors over all q-points
      integer ngqmax
! index from G+q-vectors to G-vectors
      integer, allocatable :: igqig(:,:)
! index from G-vectors to G+q-vectors
      integer, allocatable :: igigq(:,:)
! G+q-vectors in lattice coordinates
      real(8), allocatable :: vgql(:,:,:)
! G+q-vectors in Cartesian coordinates
      real(8), allocatable :: vgqc(:,:,:)
! Transformation matrix between IPW's and OIPW's
      complex(8), allocatable :: sgi(:,:)        
! Matrix element between IPW's and PW's       
      complex(8), allocatable :: mpwipw(:,:) 

!-------------------------------------------------------------------!
!     Matrix representation of the symmetry operations in MB basis  !
!-------------------------------------------------------------------!
      complex(8), allocatable :: rotmat(:,:)

!----------------------------!
!     Bare Coulomb matrix !
!----------------------------!
! G-vector cut-off for the bare coulomb matrix
      real(8) :: gmaxbarc
! number of G-vectors for the bare coulomb matrix
      integer(4), allocatable :: ngbarc(:)
! maximum number of G-vectors for the bare coulomb matrix      
      integer(4) :: ngbarcmax
! index from G+q-vectors to G-vectors for the coulomb Gmax cutoff
      integer, allocatable :: igqigb(:,:)
! index from G-vectors to G+q-vectors for the coulomb Gmax cutoff
      integer, allocatable :: igigqb(:,:)
! parameter eta used for the calculation of the structure constants      
      real(8) :: eta
! cutoff for the real space summation      
      real(8) :: rcf 
! cutoff for the reciprocal space summation           
      real(8) :: gcf 
! The lattice summations matrix      
      complex(8), allocatable :: sgm(:,:,:) 
! The matrix representation of the bare coulomb potentianl in the mixed basis            
      complex(8), allocatable :: barc(:,:) 
! The matrix representation of the square root of the bare coulomb potential
       complex(8), allocatable :: sqbarc(:,:)
! The matrix representation of the singular term of the bare coulomb matrix
      complex(8), allocatable :: barcs(:,:) 
! The square root of the singular term of the bare coulomb matrix
      complex(8), allocatable :: sqbarcs(:,:) 

! the number of eigenvectors of bare Coulomb matrix used as basis functions 
      integer(4) :: mbsiz     
! full set of the eigenvalues of barcoul matrix
      real(8), allocatable :: barcev(:)
! full set of eigenvectors of barcoul matrix        
      complex(8), allocatable:: vmat(:,:)
! reduced set of eigenvectors of barcoul matrix after barcevtol
      complex(8), allocatable:: vbas(:,:)
! transform matrix that diagonalized original bare Coulomb matrix 
      complex(8), allocatable :: barcvm(:,:)    
      
!----------------------------!
!     Dielectric matrix      !
!----------------------------!

! The dielectric function matrix matrix
      complex(8), allocatable :: epsilon(:,:,:)
! The inverse of the dielectric matrix
      complex(8), allocatable :: inveps(:,:,:) 
! Head of the dielectric function
      complex(8), allocatable :: head(:)      
! the vertical   wing of the dielectric matrix  
      complex(8), allocatable :: epsw1(:,:)
! the horizontal wing of the dielectric matrix
      complex(8), allocatable :: epsw2(:,:) 
! the macroscopic dielectric function  
      complex(8), allocatable :: emac(:,:) ! emac(1,:) -- without local field effect     
                                           ! emac(2,:) -- with    local field effect
 
      real(8) :: q0eps(3)

!---------------------------------!
!     momentum matrix elements    !
!---------------------------------!
! radial integrals coefficients (APW-core)
      real(8), allocatable :: ripacor(:,:,:,:,:,:)
! radial integrals coefficients (core-APW)
      real(8), allocatable :: ripcora(:,:,:,:,:,:)
      
!----------------------------!
! Screened coulomb potential !
!----------------------------!
      complex(8), allocatable :: wxc(:,:,:)       
! Singular term (q^-1) of the screened coulomb potential       
      complex(8), allocatable :: ws1(:,:,:)       
! Singular term (q^-2) of the screened coulomb potential       
      complex(8), allocatable :: ws2(:,:,:)     

!----------------------------!
!     Self energy            ! 
!----------------------------!
! Option for calculating selfenergy
!          0 -- perturbative G0W0 without energy shift
!          1 -- perturbative G0W0 with energy shift
!          2 -- iterative G0W0 with energy shift
!          3 -- iterative G0W0 without energy shift
      integer(4) :: iopes
! The correlation term of the selfenergy
      complex(8), allocatable :: selfec(:,:,:) 
! The exchange term of the selfenergy
      complex(8), allocatable :: selfex(:,:) 
! Singular term (q^-2) of the exchange selfenergy
      complex(8), allocatable :: selfxs2(:,:)
! Correction factors for (q^-1) and (q^-2) singularities
      real(8) :: singc1,singc2
      
!----------------------------!
!      analytic continuation !
!----------------------------!
! Type of analytic continuation
      integer(4) :: iopac ! 1/2
! Number of poles used in AC
      integer(4) :: npol
! parameters of the functions fitting the selfenergy
      complex(8), allocatable :: sacpar(:,:,:)
      
!----------------------------!
!     qp energies            !
!----------------------------!
! qp energies
      real(8), allocatable :: eqp(:,:)

!---------------------------------------!
!     Macroscopic dielectric function   !
!---------------------------------------!
! Input energies for the dielectric function (lda or gwa)
      character(3) :: epsin

!----------------------------!
!     Bandstructure          !
!----------------------------!
! GW data 
      integer(4) :: nkp1
      real(8), allocatable :: kvecs1(:,:)
      real(8), allocatable :: eks1(:,:), eqp1(:,:)
! Data from excitings bandstructure calculations
      integer(4) :: nkp2
      real(8), allocatable :: kvecs2(:,:)
      real(8), allocatable :: eks2(:,:), eqp2(:,:)

!--------------------------------------------!
!    The variables needed for angular 
!    integrals in the unit sphere
!--------------------------------------------!
      complex(8), allocatable :: sphar(:,:) !  Value of the 
                                            !  spherical harmonic
                                            !  lm (ll=l^2+l+m+1) 
                                            !  at the gridpoint 
                                            !  ileb.

      real(8), allocatable :: wleb(:) ! weight of the grid point ileb.
      integer(4) :: nleb              ! number of grip point in the sphere.

contains

      logical function gammapoint(iq)
        use mod_qpoint, only: vql
        implicit none
        integer(4) :: iq
        real(8)    :: len

        len=vql(1,iq)**2+vql(2,iq)**2+vql(3,iq)**2
        if (len.gt.1.0d-6) then
          gammapoint=.false.
        else
          gammapoint=.true.
        endif
        
        return
      end function gammapoint

end module modgw
!EOP
