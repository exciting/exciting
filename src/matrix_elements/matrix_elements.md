## Specifying the operator

Inside the muffin-tin sphere (MT) of atom \(\alpha\), the local operator must be given
either as an expansion in complex spherical harmonics
\[ O({\bf r}_\alpha) = \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^l o^{\alpha}_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
or in real spherical harmonics
\[ O({\bf r}_\alpha) = \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^l o^{\alpha}_{lm}(r_\alpha) \, R_{lm}(\hat{\bf r}_\alpha)\;, \]
where \(o^{\alpha}_{lm}(r)\) are either real or complex valued radial functions.

Inside the IR, the operator must be given as a Fourier series
\[ O({\bf r}) = \sum_{\bf G+q} \hat{O}({\bf G+q}) \, {\rm e}^{{\rm i}({\bf G+q}) \cdot {\bf r}} \;.\]
Note, that the operator can be a Bloch wave with wavevector \({\bf q}\) and must be provided
on the real-space FFT grid corresponding to the set of \({\bf G}\)-vectors it is expanded in.
If \({\bf q} \neq 0\), only the cell-periodic part of the operator must be provided.

## Usage of the module

The calculation of matrix elements is made up of the following steps:

1.  Initialization of the module using the subroutine [[me_init(subroutine)]].
2.  Preparation of the MT and IR parts that are independent of the \({\bf k}\)-vector
    of the basis functions / wavefunctions using the subroutines [[me_mt_prepare(subroutine)]]
    and [[me_ir_prepare(subroutine)]].
3.  Loop over \({\bf k}\)-vectors and actual evaluation of the matrix elements 
    using the subroutines [[me_mt_mat(subroutine)]] and [[me_ir_mat(subroutine)]].
4.  Finalization of the module using the subroutine [[me_finit(subroutine)]].

In the following subsection, a code example for the usage of this module is presented.

### Code example

This examples shows how to compute the matrix elements of the local effective potential
which is, besides the kinetic energy term, the fundamental contribution to the Kohn-Sham Hamiltonian.

The following variables are assumed to be provided and initialized in the scope of the usage of the module
(names of code globals are used but actual names and values can be arbitray):

```fortran
! set of k-vectors
type(k_set) :: kset
! set of G-vectors the potential is expanded in
type(G_set) :: Gset
! set of G+k-vectors for (L)APW basis functions (generated using `Gset` and `kset`)
type(Gk_set) :: Gkset
! maximum l for potential expansion
integer :: lmaxvr
! `mt_basis_type` object carrying the (L)APW and LO radial functions
type(mt_basis_type) :: mt_basis
! effective potential in MTs given in a real spherical harmonics expansion
real(dp), allocatable :: veffmt(:,:,:)
! effective potential in IR given on the real space FFT grid corresponding to `Gset`
real(dp), allocatable :: veffir(:)
```

The first step is the declaration of additional variables and the preparation of the module.
Here, we also employ the subroutines [[me_mt_alloc(subroutine)]] and [[me_ir_alloc(subroutine)]].

```fortran
! Hamiltonian matrix
complex(dp), allocatable :: H(:,:)
! radial integrals times Gaunt coefficients
! Note: `rigntHmt` is complex because the potential is given by
! real spherical harmonics and as a consequence, the Gaunt coefficients
! of type <Y|R|Y> are complex.
complex(dp), allocatable :: rigntHmt(:,:,:)
! interstitial potential times characteristic function in reciprocal space
complex(dp), allocatable :: veffig(:)
! (L)APW matching coefficients
complex(dp), allocatable :: apwalm(:,:,:,:)

! allocate Hamiltonian matrix to maximum matrix size
allocate( H(Gkset%ngkmax+nlotot, Gkset%ngkmax+nlotot) )

! allocate (L)APW matching coefficients
allocate( apwalm(Gkset%ngkmax, apwordmax, lmmaxapw, natmtot) )

! initialize matrix_elements module
! Note: Must be called before `me_XX_alloc`.
call me_init( mt_basis, lmaxvr, Gset )

! allocate radial integrals times Gaunt coefficients
call me_mt_alloc( rigntHmt )

! allocate interstitial potential times characteristic function in reciprocal space
call me_ir_alloc( veffig )
```

In the second step, the \({\bf k}\)-point independent variables are initialized which are
the radial integrals times the Gaunt coefficients inside the MTs and the interstitial
operator times the characteristic function in reciprocal space.

```fortran
! calculate radial MT integrals times Gaunt coefficients for each atom
do is = 1, nspecies
  do ia = 1, natoms(is)
    ias = idxas(ia, is)

    call me_mt_prepare( is, ia, lmaxvr, zone, veffmt(:, :, ias), zzero, rigntHmt(:, :, ias) )
    
  end do
end do

! calculate interstitial potential times characteristic function in reciprocal space
call me_ir_prepare( zone, veffir, zzero, veffig )
```

In the third step, for each \({\bf k}\)-point, the matrix elements are calculated.

```fortran
do ik = 1, kset%nkpt
  ! set matrix elements to zero
  ! Note: The contributions from each muffin-tin and the interstitial
  !       will be added to the input argument `H`.
  H = zzero

  ! get (L)APW matching coefficients
  call match( Gkset%ngk(1, ik), Gkset%gkc(:, 1, ik), Gkset%tpgkc(:, :, 1, ik), Gkset%sfacgk(:, :, 1, ik), &
              apwalm )

  ! add contribution coming from each MT
  do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia, is)
      call me_mt_mat( is, ia, Gkset%ngk(1, ik), apwalm(:, :, :, ias), zone, rigntHmt(:, :, ias), zone, H )
    end do
  end do

  ! add contribution coming from IR
  call me_ir_mat( Gkset, ik, zone, veffig, zone, H )

  ! Here, additional contributions to the Hamiltonian might be added (e.g. kinetic energy)
  ! and diagonalization might follow.
end do
```

At last, after the calculation of matrix elements is finished, memory can be freed
from unneeded variables.

```fortran
deallocate( H, rigntHmt, veffig, apwalm )
call me_finit
```

See also [[matrix_elements_test(module)]] for further application.
There, the module is used to compute matrix elements of plane waves.

## Muffin-tin representation of basis functions

For a detailed description of the (L)APW+LO basis representation inside the muffin-tin spheres,
see [[muffin_tin_basis(module)]].
