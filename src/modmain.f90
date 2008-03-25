
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmain

!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
real(8) avec(3,3)
! inverse of lattice vector matrix
real(8) ainv(3,3)
! reciprocal lattice vectors
real(8) bvec(3,3)
! inverse of reciprocal lattice vector matrix
real(8) binv(3,3)
! unit cell volume
real(8) omega
! any vector with length less than epslat is considered zero
real(8) epslat

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=200
! number of species
integer nspecies
! number of atoms for each species
integer natoms(maxspecies)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot
! index to atoms and species
integer idxas(maxatoms,maxspecies)
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! atomic positions in lattice coordinates
real(8) atposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc(3,maxatoms,maxspecies)

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species file names
character(256) spfname(maxspecies)
! species name
character(256) spname(maxspecies)
! species symbol
character(256) spsymb(maxspecies)
! species nuclear charge
real(8) spzn(maxspecies)
! species electronic charge
real(8) spze(maxspecies)
! species mass
real(8) spmass(maxspecies)
! smallest radial point for each species
real(8) sprmin(maxspecies)
! effective infinity for species
real(8) sprmax(maxspecies)
! number of radial points to effective infinity for each species
integer spnr(maxspecies)
! maximum spnr over all the species
integer spnrmax
! maximum allowed states for each species
integer, parameter :: maxspst=40
! number of states for each species
integer spnst(maxspecies)
! maximum spnst over all the species
integer spnstmax
! state principle quantum number for each species
integer spn(maxspst,maxspecies)
! state l value for each species
integer spl(maxspst,maxspecies)
! state k value for each species
integer spk(maxspst,maxspecies)
! spcore is .true. if species state is core
logical spcore(maxspst,maxspecies)
! state eigenvalue for each species
real(8) speval(maxspst,maxspecies)
! state occupancy for each species
real(8) spocc(maxspst,maxspecies)
! species radial mesh
real(8), allocatable :: spr(:,:)
! species charge density
real(8), allocatable :: sprho(:,:)
! species self-consistent potential
real(8), allocatable :: spvr(:,:)

!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! radial function integration and differentiation polynomial order
integer nprad
! number of muffin-tin radial points for each species
integer nrmt(maxspecies)
! maximum nrmt over all the species
integer nrmtmax
! autormt is .true. for automatic determination of muffin-tin radii
logical autormt
! parameters for determining muffin-tin radii automatically
real(8) rmtapm(2)
! muffin-tin radii
real(8) rmt(maxspecies)
! smallest muffin-tin radius
real(8) rmtmin
! largest muffin-tin radius
real(8) rmtmax
! radial step length for coarse mesh
integer lradstp
! number of coarse radial mesh points
integer nrcmt(maxspecies)
! maximum nrcmt over all the species
integer nrcmtmax
! coarse muffin-tin radial mesh
real(8), allocatable :: rcmt(:,:)
! maximum allowable angular momentum for augmented plane waves
integer, parameter :: maxlapw=50
! maximum angular momentum for augmented plane waves
integer lmaxapw
! (lmaxapw+1)^2
integer lmmaxapw
! maximum angular momentum for potentials and densities
integer lmaxvr
! (lmaxvr+1)^2
integer lmmaxvr
! maximum angular momentum used when evaluating the Hamiltonian matrix elements
integer lmaxmat
! (lmaxmat+1)^2
integer lmmaxmat
! fraction of muffin-tin radius which constitutes the inner part
real(8) fracinr
! maximum angular momentum in the inner part of the muffin-int
integer lmaxinr
! (lmaxinr+1)^2
integer lmmaxinr
! number of radial points to the inner part of the muffin-tin
integer nrmtinr(maxspecies)
! index to (l,m) pairs
integer, allocatable :: idxlm(:,:)

!--------------------------------!
!     spin related variables     !
!--------------------------------!
! spinpol is .true. for spin-polarised calculations
logical spinpol
! spinorb is .true. for spin-orbit coupling
logical spinorb
! fixspin type: 0 = none, 1 = global, 2 = local, 3 = global + local
integer fixspin
! dimension of magnetisation and magnetic vector fields (1 or 3)
integer ndmag
! fixed total spin magnetic moment
real(8) momfix(3)
! fixed spin moment global effective field in Cartesian coordinates
real(8) bfsmc(3)
! muffin-tin fixed spin moments
real(8) mommtfix(3,maxatoms,maxspecies)
! muffin-tin fixed spin moment effective fields in Cartesian coordinates
real(8) bfsmcmt(3,maxatoms,maxspecies)
! fixed spin moment field mixing parameter
real(8) taufsm
! second-variational spinor dimension (1 or 2)
integer nspinor
! external magnetic field in each muffin-tin in lattice coordinates
real(8) bflmt(3,maxatoms,maxspecies)
! external magnetic field in each muffin-tin in Cartesian coordinates
real(8) bfcmt(3,maxatoms,maxspecies)
! global external magnetic field in lattice coordinates
real(8) bfieldl(3)
! global external magnetic field in Cartesian coordinates
real(8) bfieldc(3)
! external magnetic fields are multiplied by reducebf after each iteration
real(8) reducebf
! spinsprl if .true. if a spin-spiral is to be calculated
logical spinsprl
! number of spin-dependent first-variational functions per state
integer nspnfv
! spin-spiral q-vector in lattice coordinates
real(8) vqlss(3)
! spin-spiral q-vector in Cartesian coordinates
real(8) vqcss(3)

!----------------------------!
!     symmetry variables     !
!----------------------------!
! nosym is .true. if no symmetry information should be used
logical nosym
! number of Bravais lattice point group symmetries
integer nsymlat
! Bravais lattice point group symmetries
integer symlat(3,3,48)
! determinants of lattice symmetry matrices (1 or -1)
integer symlatd(48)
! index to inverses of the lattice symmetries
integer isymlat(48)
! lattice point group symmetries in Cartesian coordinates
real(8) symlatc(3,3,48)
! tshift is .true. if atomic basis is allowed to be shifted
logical tshift
! maximum of symmetries allowed
integer, parameter :: maxsymcrys=192
! number of crystal symmetries
integer nsymcrys
! crystal symmetry translation vector in lattice coordinates
real(8) vtlsymc(3,maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
integer lsplsymc(maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
integer lspnsymc(maxsymcrys)
! equivalent atom index for each crystal symmetry
integer, allocatable :: ieqatom(:,:,:)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
logical, allocatable :: eqatoms(:,:,:)
! number of site symmetries
integer, allocatable :: nsymsite(:)
! site symmetry spatial rotation element in lattice point group
integer, allocatable :: lsplsyms(:,:)
! site symmetry global spin rotation element in lattice point group
integer, allocatable :: lspnsyms(:,:)

!--------------------------------!
!     G-vector set variables     !
!--------------------------------!
! G-vector cut-off for interstitial potential and density
real(8) gmaxvr
! G-vector grid sizes
integer ngrid(3)
! total number of G-vectors
integer ngrtot
! integer grid intervals for each direction
integer intgv(3,2)
! number of G-vectors with G < gmaxvr
integer ngvec
! locations of G-vectors on integer grid
integer, allocatable :: ivg(:,:)
! map from integer grid to G-vector array
integer, allocatable :: ivgig(:,:,:)
! map from G-vector array to FFT array
integer, allocatable :: igfft(:)
! G-vectors in Cartesian coordinates
real(8), allocatable :: vgc(:,:)
! length of G-vectors
real(8), allocatable :: gc(:)
! spherical harmonics of the G-vectors
complex(8), allocatable :: ylmg(:,:)
! structure factor for the G-vectors
complex(8), allocatable :: sfacg(:,:)
! G-space characteristic function: 0 inside the muffin-tins and 1 outside
complex(8), allocatable :: cfunig(:)
! real-space characteristic function: 0 inside the muffin-tins and 1 outside
real(8), allocatable :: cfunir(:)
! damping coefficient for characteristic function
real(8) cfdamp

!-------------------------------!
!     k-point set variables     !
!-------------------------------!
! maximum de Broglie wavelength
real(8) rlambda
! autokpt is .true. if the k-point set is determined by rlambda
logical autokpt
! k-point grid sizes
integer ngridk(3)
! total number of k-points
integer nkpt
! k-point offset
real(8) vkloff(3)
! reducek is .true. if k-points are to be reduced (with crystal symmetries)
logical reducek
! locations of k-points on integer grid
integer, allocatable :: ivk(:,:)
! k-points in lattice coordinates
real(8), allocatable :: vkl(:,:)
! k-points in Cartesian coordinates
real(8), allocatable :: vkc(:,:)
! k-point weights
real(8), allocatable :: wkpt(:)
! map from non-reduced grid to reduced set
integer, allocatable :: ikmap(:,:,:)
! total number of non-reduced k-points
integer nkptnr
! locations of non-reduced k-points on integer grid
integer, allocatable :: ivknr(:,:)
! non-reduced k-points in lattice coordinates
real(8), allocatable :: vklnr(:,:)
! non-reduced k-points in Cartesian coordinates
real(8), allocatable :: vkcnr(:,:)
! non-reduced k-point weights
real(8), allocatable :: wkptnr(:)
! map from non-reduced grid to non-reduced set
integer, allocatable :: ikmapnr(:,:,:)
! k-point at which to determine effective mass tensor
real(8) vklem(3)
! displacement size for computing the effective mass tensor
real(8) deltaem
! number of displacements in each direction
integer ndspem

!----------------------------------!
!     G+k-vector set variables     !
!----------------------------------!
! smallest muffin-tin radius times gkmax
real(8) rgkmax
! maximum |G+k| cut-off for APW functions
real(8) gkmax
! number of G+k-vectors for augmented plane waves
integer, allocatable :: ngk(:,:)
! maximum number of G+k-vectors over all k-points
integer ngkmax
! index from G+k-vectors to G-vectors
integer, allocatable :: igkig(:,:,:)
! G+k-vectors in lattice coordinates
real(8), allocatable :: vgkl(:,:,:,:)
! G+k-vectors in Cartesian coordinates
real(8), allocatable :: vgkc(:,:,:,:)
! length of G+k-vectors
real(8), allocatable :: gkc(:,:,:)
! (theta, phi) coordinates of G+k-vectors
real(8), allocatable :: tpgkc(:,:,:,:)
! structure factor for the G+k-vectors
complex(8), allocatable :: sfacgk(:,:,:,:)

!-------------------------------!
!     q-point set variables     !
!-------------------------------!
! q-point grid sizes
integer ngridq(3)
! total number of q-points
integer nqpt
! reduceq is .true. if q-points are to be reduced (with crystal symmetries)
logical reduceq
! locations of q-points on integer grid
integer, allocatable :: ivq(:,:)
! map from non-reduced grid to reduced set
integer, allocatable :: iqmap(:,:,:)
! q-points in lattice coordinates
real(8), allocatable :: vql(:,:)
! q-points in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! q-point weights
real(8), allocatable :: wqpt(:)
! weights associated with the integral of 1/q^2
real(8), allocatable :: wiq2(:)

!-----------------------------------------------------!
!     spherical harmonic transform (SHT) matrices     !
!-----------------------------------------------------!
! real backward SHT matrix for lmaxapw
real(8), allocatable :: rbshtapw(:,:)
! real forward SHT matrix for lmaxvr
real(8), allocatable :: rfshtvr(:,:)
! complex backward SHT matrix for lmaxapw
complex(8), allocatable :: zbshtapw(:,:)
! complex forward SHT matrix for lmaxapw
complex(8), allocatable :: zfshtapw(:,:)
! complex forward SHT matrix for lmaxvr
complex(8), allocatable :: zfshtvr(:,:)

!-----------------------------------------!
!     potential and density variables     !
!-----------------------------------------!
! exchange-correlation functional type
integer xctype
! exchange-correlation functional description
character(256) xcdescr
! exchange-correlation functional spin treatment
integer xcspin
! exchange-correlation functional density gradient treatment
integer xcgrad
! muffin-tin charge density
real(8), allocatable :: rhomt(:,:,:)
! interstitial real-space charge density
real(8), allocatable :: rhoir(:)
! muffin-tin magnetisation vector field
real(8), allocatable :: magmt(:,:,:,:)
! interstitial magnetisation vector field
real(8), allocatable :: magir(:,:)
! muffin-tin Coulomb potential
real(8), allocatable :: vclmt(:,:,:)
! interstitial real-space Coulomb potential
real(8), allocatable :: vclir(:)
! order of polynomial for pseudocharge density
integer npsden
! muffin-tin exchange-correlation potential
real(8), allocatable :: vxcmt(:,:,:)
! interstitial real-space exchange-correlation potential
real(8), allocatable :: vxcir(:)
! muffin-tin exchange-correlation magnetic field
real(8), allocatable :: bxcmt(:,:,:,:)
! interstitial exchange-correlation magnetic field
real(8), allocatable :: bxcir(:,:)
! nosource is .true. if the field is to be made source-free
logical nosource
! muffin-tin effective potential
real(8), allocatable :: veffmt(:,:,:)
! interstitial effective potential
real(8), allocatable :: veffir(:)
! G-space interstitial effective potential
complex(8), allocatable :: veffig(:)
! muffin-tin exchange energy density
real(8), allocatable :: exmt(:,:,:)
! interstitial real-space exchange energy density
real(8), allocatable :: exir(:)
! muffin-tin correlation energy density
real(8), allocatable :: ecmt(:,:,:)
! interstitial real-space correlation energy density
real(8), allocatable :: ecir(:)
! default potential mixing parameter
real(8) beta0
! maximum potential mixing parameter
real(8) betamax

!-------------------------------------!
!     charge and moment variables     !
!-------------------------------------!
! tolerance for error in total charge
real(8) epschg
! total nuclear charge
real(8) chgzn
! total core charge
real(8) chgcr
! core leakage charge
real(8) chgcrlk
! total valence charge
real(8) chgval
! excess charge
real(8) chgexs
! total charge
real(8) chgtot
! calculated total charge
real(8) chgcalc
! interstitial region charge
real(8) chgir
! muffin-tin charges
real(8), allocatable :: chgmt(:)
! total muffin-tin charge
real(8) chgmttot
! total moment
real(8) momtot(3)
! interstitial region moment
real(8) momir(3)
! muffin-tin moments
real(8), allocatable :: mommt(:,:)
! total muffin-tin moment
real(8) mommttot(3)

!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! maximum allowable APW order
integer, parameter :: maxapword=3
! APW order
integer apword(0:maxlapw,maxspecies)
! maximum of apword over all angular momenta and species
integer apwordmax
! APW initial linearisation energies
real(8) apwe0(maxapword,0:maxlapw,maxspecies)
! APW linearisation energies
real(8), allocatable :: apwe(:,:,:)
! APW derivative order
integer apwdm(maxapword,0:maxlapw,maxspecies)
! apwve is .true. if the linearisation energies are allowed to vary
logical apwve(maxapword,0:maxlapw,maxspecies)
! APW radial functions
real(8), allocatable :: apwfr(:,:,:,:,:)
! derivate of radial functions at the muffin-tin surface
real(8), allocatable :: apwdfr(:,:,:)
! maximum number of local-orbitals
integer, parameter :: maxlorb=20
! maximum allowable local-orbital order
integer, parameter :: maxlorbord=4
! number of local-orbitals
integer nlorb(maxspecies)
! maximum nlorb over all species
integer nlomax
! total number of local-orbitals
integer nlotot
! local-orbital order
integer lorbord(maxlorb,maxspecies)
! local-orbital angular momentum
integer lorbl(maxlorb,maxspecies)
! maximum lorbl over all species
integer lolmax
! (lolmax+1)^2
integer lolmmax
! local-orbital initial energies
real(8) lorbe0(maxlorbord,maxlorb,maxspecies)
! local-orbital energies
real(8), allocatable :: lorbe(:,:,:)
! local-orbital derivative order
integer lorbdm(maxlorbord,maxlorb,maxspecies)
! lorbve is .true. if the linearisation energies are allowed to vary
logical lorbve(maxlorbord,maxlorb,maxspecies)
! local-orbital radial functions
real(8), allocatable :: lofr(:,:,:,:)
! energy step size for locating the band energy
real(8) deband

!-------------------------------------------!
!     overlap and Hamiltonian variables     !
!-------------------------------------------!
! order of overlap and Hamiltonian matrices for each k-point
integer, allocatable :: nmat(:,:)
! maximum nmat over all k-points
integer nmatmax
! size of packed matrices
integer, allocatable :: npmat(:,:)
! index to the position of the local-orbitals in the H and O matrices
integer, allocatable :: idxlo(:,:,:)
! APW-local-orbital overlap integrals
real(8), allocatable :: oalo(:,:,:)
! local-orbital-local-orbital overlap integrals
real(8), allocatable :: ololo(:,:,:)
! APW-APW Hamiltonian integrals
real(8), allocatable :: haa(:,:,:,:,:,:)
! local-orbital-APW Hamiltonian integrals
real(8), allocatable :: hloa(:,:,:,:,:)
! local-orbital-local-orbital Hamiltonian integrals
real(8), allocatable :: hlolo(:,:,:,:)
! complex Gaunt coefficient array
complex(8), allocatable :: gntyry(:,:,:)

!--------------------------------------------!
!     eigenvalue and occupancy variables     !
!--------------------------------------------!
! number of empty states
integer nempty
! number of first-variational states
integer nstfv
! number of second-variational states
integer nstsv
! smearing type
integer stype
! smearing function description
character(256) sdescr
! smearing width
real(8) swidth
! maximum allowed occupancy (1 or 2)
real(8) occmax
! convergence tolerance for occupancies
real(8) epsocc
! second-variational occupation number array
real(8), allocatable :: occsv(:,:)
! spin characters of the second-variational states
real(8), allocatable :: spnchr(:,:,:)
! Fermi energy for second-variational states
real(8) efermi
! density of states at the Fermi energy
real(8) fermidos
! minimum allowed eigenvalue
real(8) evalmin
! second-variational eigenvalues
real(8), allocatable :: evalsv(:,:)
! tevecsv is .true. if second-variational eigenvectors are calculated
logical tevecsv
! maximum number of k-point and states indices in user-defined list
integer, parameter :: maxkst=20
! number of k-point and states indices in user-defined list
integer nkstlist
! user-defined list of k-point and state indices
integer kstlist(2,maxkst)

!------------------------------!
!     core state variables     !
!------------------------------!
! eigenvalues for core states
real(8), allocatable :: evalcr(:,:)
! radial wavefunctions for core states
real(8), allocatable :: rwfcr(:,:,:,:)
! radial charge density for core states
real(8), allocatable :: rhocr(:,:)

!--------------------------!
!     energy variables     !
!--------------------------!
! eigenvalue sum
real(8) evalsum
! electron kinetic energy
real(8) engykn
! core electron kinetic energy
real(8) engykncr
! nuclear-nuclear energy
real(8) engynn
! electron-nuclear energy
real(8) engyen
! Hartree energy
real(8) engyhar
! Coulomb energy (E_nn + E_en + E_H)
real(8) engycl
! electronic Coulomb potential energy
real(8) engyvcl
! Madelung term
real(8) engymad
! exchange-correlation potential energy
real(8) engyvxc
! exchange-correlation effective field energy
real(8) engybxc
! energy of external global magnetic field
real(8) engybext
! energy of muffin-tin magnetic fields (non-physical)
real(8) engybmt
! exchange energy
real(8) engyx
! correlation energy
real(8) engyc
! compensating background charge energy
real(8) engycbc
! total energy
real(8) engytot

!-------------------------!
!     force variables     !
!-------------------------!
! tforce is .true. if force should be calculated
logical tforce
! tfibs is .true. if the IBS contribution to the force is to be calculated
logical tfibs
! Hellmann-Feynman force on each atom
real(8), allocatable :: forcehf(:,:)
! core correction to force on each atom
real(8), allocatable :: forcecr(:,:)
! IBS core force on each atom
real(8), allocatable :: forceibs(:,:)
! total force on each atom
real(8), allocatable :: forcetot(:,:)
! previous total force on each atom
real(8), allocatable :: forcetp(:,:)
! maximum force magnitude over all atoms
real(8) forcemax
! default step size parameter for structural optimisation
real(8) tau0atm
! step size parameters for each atom
real(8), allocatable :: tauatm(:)

!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
integer maxscl
! current self-consistent loop number
integer iscl
! effective potential convergence tolerance
real(8) epspot
! energy convergence tolerance
real(8) epsengy
! force convergence tolerance
real(8) epsforce

!------------------------------------------------!
!     density of states and optics variables     !
!------------------------------------------------!
! number of energy intervals in the DOS/optics function
integer nwdos
! effective size of k/q-point grid for integrating the Brillouin zone
integer ngrdos
! smoothing level for DOS/optics function
integer nsmdos
! energy interval for DOS/optics function
real(8) wdos(2)
! scissors correction
real(8) scissor
! number of optical matrix components required
integer noptcomp
! required optical matrix components
integer optcomp(3,27)
! usegdft is .true. if the generalised DFT correction is to be used
logical usegdft
! intraband is .true. if the intraband term is to be added to the optical matrix
logical intraband
! bcsym is .true. if the band characters are to correspond to the the
! irreducible representations of the site symmetries
logical bcsym

!-------------------------------------!
!     1D/2D/3D plotting variables     !
!-------------------------------------!
! number of vertices in 1D plot
integer nvp1d
! total number of points in 1D plot
integer npp1d
! vertices in lattice coordinates for 1D plot
real(8), allocatable :: vvlp1d(:,:)
! distance to vertices in 1D plot
real(8), allocatable :: dvp1d(:)
! plot vectors in lattice coordinates for 1D plot
real(8), allocatable :: vplp1d(:,:)
! distance to points in 1D plot
real(8), allocatable :: dpp1d(:)
! corner vectors of 2D plot in lattice coordinates
real(8) vclp2d(3,3)
! grid sizes of 2D plot
integer np2d(2)
! number of units cells in 3D plot
integer nup3d(3)
! grid sizes in 3D plot
integer np3d(3)
! number of states for plotting Fermi surface
integer nstfsp

!-------------------------------------------------------------------!
!     non-local matrix elements, OEP and Hartree-Fock variables     !
!-------------------------------------------------------------------!
! maximum number of core states over all species
integer ncrmax
! maximum number of OEP iterations
integer maxitoep
! default and maximum step size and mixing parameter for OEP
real(8) tauoep(3)
! magnitude of the OEP residue
real(8) resoep
! kinetic matrix elements
complex(8), allocatable :: kinmatc(:,:,:)
! complex versions of the exchange potential and field
complex(8), allocatable :: zvxmt(:,:,:)
complex(8), allocatable :: zvxir(:)
complex(8), allocatable :: zbxmt(:,:,:,:)
complex(8), allocatable :: zbxir(:,:)

!-------------------------!
!     LDA+U variables     !
!-------------------------!
! type of LDA+U to use (0: none)
integer ldapu
! maximum angular momentum
integer, parameter :: lmaxlu=3
integer, parameter :: lmmaxlu=(lmaxlu+1)**2
! angular momentum for each species
integer llu(maxspecies)
! U and J values for each species
real(8) ujlu(2,maxspecies)
! LDA+U density matrix
complex(8), allocatable :: dmatlu(:,:,:,:,:)
! LDA+U potential matrix in (l,m) basis
complex(8), allocatable :: vmatlu(:,:,:,:,:)
! LDA+U energy for each atom
real(8), allocatable :: engyalu(:)
! interpolation constant alpha for each atom (PRB 67, 153106 (2003))
real(8), allocatable :: alphalu(:)
! energy from the LDA+U correction
real(8) engylu

!--------------------------!
!     phonon variables     !
!--------------------------!
! number of primitive unit cells in phonon supercell
integer nphcell
! phonon displacement distance
real(8) deltaph
! number of vectors for writing out frequencies and eigenvectors
integer nphwrt
! vectors in lattice coordinates for writing out frequencies and eigenvectors
real(8), allocatable :: vqlwrt(:,:)

!-------------------------------------------------------------!
!     reduced density matrix functional (RDMFT) variables     !
!-------------------------------------------------------------!
! non-local matrix elements for varying occupation numbers
real(8), allocatable :: vnlrdm(:,:,:,:)
! Coulomb potential matrix elements
complex(8), allocatable :: vclmat(:,:,:)
! derivative of kinetic energy w.r.t. natural orbital coefficients
complex(8), allocatable :: dkdc(:,:,:)
! step size for occupation numbers
real(8) taurdmn
! step size for natural orbital coefficients
real(8) taurdmc
! xc functional
integer rdmxctype
! maximum number of self-consistent loops
integer rdmmaxscl
! maximum number of iterations for occupation number optimisation
integer maxitn
! maximum number of iteration for natural orbital optimisation
integer maxitc
! exponent for the functional
real(8) rdmalpha

!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
real(8) timeinit
! Hamiltonian and overlap matrix set up
real(8) timemat
! first-variational calculation
real(8) timefv
! second-variational calculation
real(8) timesv
! charge density calculation
real(8) timerho
! potential calculation
real(8) timepot
! force calculation
real(8) timefor

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! square root of two
real(8), parameter :: sqtwo=1.4142135623730950488d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
! complex constants
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zhalf=(0.5d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0)
! array of i**l values
complex(8), allocatable :: zil(:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
complex(8) sigmat(2,2,3)
data sigmat / (0.d0,0.d0), (1.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
              (0.d0,0.d0), (0.d0,1.d0),(0.d0,-1.d0), (0.d0,0.d0), &
              (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0),(-1.d0,0.d0) /

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer version(3)
data version / 0,9,151 /
! maximum number of tasks
integer, parameter :: maxtasks=20
! number of tasks
integer ntasks
! task array
integer tasks(maxtasks)
! current task
integer task
! tstop is .true. if STOP file exists
logical tstop
! tlast is .true. if self-consistent loop is on the last iteration
logical tlast
! number of iterations after which STATE.OUT is written
integer nwrite
! file name extension for files generated by gndstate
character(256) filext
! default file extension
data filext / '.OUT' /
! scratch space path
character(256) scrpath
! maximum number of note lines
integer, parameter :: maxnlns=20
! number of note lines
integer notelns
! notes to include in INFO.OUT
character(80) notes(maxnlns)

end module

