!> title:   muffin tin basis
!> summary: Module for generalization and convenient use of 
!>          mixed (L)APW and LO basis functions inside the muffin tins.
!> author:  Sebastian Tillack
!> date:    January 2023
!> licence: GPL
!>
!> Inside the muffin-tin sphere of atom \(\alpha\) the (L)APW basis functions are given by
!> \[ \phi_{\bf G+k}({\bf r}) = \sum_{l,m,\xi} A^{\alpha}_{lm,\xi,{\bf G+k}} \, 
!> u^\alpha_{l,\xi}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha)\;, \]
!> where \(\xi\) is the APW linearization order, \(A^{\alpha}_{lm,\xi,{\bf G+k}}\) are the
!> APW matching coefficients, \({\bf r}_\alpha = {\bf r} - {\bf \tau}_\alpha \) with 
!> \({\bf \tau}_\alpha\) being the atom position, \(u^\alpha_{l,\xi}(r_\alpha)\) are the radial
!> functions and \(Y_{lm}(\hat{\bf r}_\alpha)\) are spherical harmonics.
!> 
!> The local orbitals are given by
!> \[ \phi_L({\bf r}) = \delta_{\alpha,\alpha_L} \, f^L(r_\alpha) \, Y_{l_L m_L}(\hat{\bf r}_\alpha) \;, \]
!> where \(L\) is the index labeling the LO, \(\alpha_L\) labels the atom, the LO belongs to
!> and \(l_L\) and \(m_L\) label its angular and magnetic quantum number, respecively.
!> 
!> A wavefunction \(\psi_{n{\rm k}}({\bf r})\) inside the muffin-tin \(\alpha\) is given by
!> \[ \psi^\alpha_{n{\bf k}}({\bf r}) = \sum_{\bf G+k} C^{n{\bf k}}_{\bf G+k}\, \phi_{\bf G+k}({\bf r}) 
!>    + \sum_L C^{n{\bf k}}_L\, \phi_L({\bf r}) \; ,\]
!> where \(C^{n{\bf k}}\) is the corresponding eigenvector of the Hamiltonian. We can rewrite the same 
!> wavefunction as 
!> \[ \psi^\alpha_{n{\bf k}}({\bf r}) = \sum_\lambda C^{n{\bf k}\alpha}_\lambda\, \phi^\alpha_\lambda({\bf r}) = 
!>    \sum_\lambda C^{n{\bf k}\alpha}_\lambda\, g^\alpha_\lambda(r_\alpha)\, Y_{l_\lambda m_\lambda}(\hat{\bf r}_\alpha)\;, \]
!> where \(\lambda\) counts all muffin-tin basis functions 
!> \({\phi^\alpha_\lambda({\bf r}) = g^\alpha_\lambda(r_\alpha)\, Y_{l_\lambda m_\lambda}(\hat{\bf r}_\alpha)}\) 
!> and is a combined \((l,m,\xi)\) index for (L)APWs and \(L\) for all LOs. The respective eigenvectors are then given by
!> \[ C^{n{\bf k}\alpha}_\lambda = \begin{cases} 
!>    \sum_{\bf G+k} C^{n{\bf k}}_{\bf G+k}\, A^\alpha_{lm, \xi, {\bf G+k}} & \lambda \rightarrow \text{ (L)APW} \\
!>    C^{n{\bf k}}_L\, \delta_{\alpha \alpha_L} & \lambda \rightarrow \text{ LO} \; . \end{cases} \]
!> The radial functions \(g^\alpha_\lambda(r)\), accordingly, are either (L)APW or LO radial functions
!> \[ g^\alpha_\lambda(r) = \begin{cases} u^\alpha_{l,\xi}(r) & \lambda \rightarrow \text{ (L)APW} \\
!>    f^L(r) & \lambda \rightarrow \text{ LO} \; . \end{cases} \]
!> Note that the radial functions \(g^\alpha_\lambda\) only depend on \(l_\lambda\) but not on \(m_\lambda\).
!> Hence, we introduce the index \(\tilde{\lambda}\) which runs over all (L)APWs and LOs with the same
!> angular momentum \(l\) and write \(g^\alpha_\tilde{\lambda}\) instead.
module muffin_tin_basis
  use precision, only: dp

  implicit none
  private

  !> Object describing the muffin-tin basis functions
  !>
  type :: mt_basis_type
    private
    !> (L)APW radial functions \(u^\alpha_{l,\xi}(r)\),
    !> The index order is `apw_rad_fun(ir, d, o, l, ias)` with radial grid point `ir`, radial derivative order `d`, 
    !> linearisation order \(\xi=\)`o`, angular momentum \(l=\)`l` and atom index \(\alpha=\)`ias`.
    real(dp), allocatable :: apw_rad_fun(:,:,:,:,:)
    !> LO radial functions \(f^L(r)\),
    !> The index order is `lo_rad_fun(ir, d, ilo, ias)` with radial grid point `ir`, radial derivative order `d`, 
    !> local orbital index \(L=\)`ilo` and atom index \(\alpha=\)`ias`.
    real(dp), allocatable :: lo_rad_fun(:,:,:,:)
    !> radial grids on which the radial functions are given for each species
    real(dp), public, allocatable :: rad_grid(:,:)
    !> number of radial grid points for each species
    integer, public, allocatable :: n_rad_grid(:)
    !> For a given species \(\alpha=\)`is`, `n_rad_fun(l, is)` gives the total number of radial functions \(g^\alpha_\tilde{\lambda}(r)\)
    !> with angular momentum \(l=\)`l`.
    integer, public, allocatable :: n_rad_fun(:,:)
    !> maximum `n_rad_fun` over all angular momenta and species
    integer, public :: n_rad_fun_max
    !> For a given species \(\alpha=\)`is`, `n_basis_fun(is)` gives the total number of basis functions \(\phi^\alpha_{\lambda}({\bf r})\).
    integer, public, allocatable :: n_basis_fun(:)
    !> maximum `n_basis_fun` over all species
    integer, public :: n_basis_fun_max
    !> maximum angular momentum \(l\) considered for basis functions \(\phi^\alpha_\lambda({\bf r})\)
    integer, public :: lmax_basis
    !> For a given combined index \((l,m)=\)`lm`, radial function index \(\tilde{\lambda}=\)`lam` and species \(\alpha=\)`is`, 
    !> `idx_basis_fun(lm, lam, is)` gives the index \(\lambda\) of the corresponding basis function \(\phi^\alpha_{\lambda}({\bf r})\).
    integer, public, allocatable :: idx_basis_fun(:,:,:)
    !> For a given radial function index \(\tilde{\lambda}=\)`lam`, angular momentum \(l=\)`l` and species \(\alpha=\)`is`, 
    !> `rad_fun_to_apw_lo(lam, l, is)` either gives the (L)APW linearization order \(\xi\) when the value is positive,
    !> or the corresponding (negative) LO index `ilo` when the value is negative. 
    !> (Intended for internal use olny. Use [[get_rad_fun(subroutine)]] to get the corresponding radial function directly.)
    integer, public, allocatable :: rad_fun_to_apw_lo(:,:,:)

    contains
      private
      !> Free memory from `mt_basis_type` type object.
      procedure, public :: destroy
      !> Get the radial function \(g^\alpha_\tilde{\lambda}(r)\) for a given index \(\tilde{\lambda}=\)`lam`, 
      !> atom index \(\alpha=\)`ias` and angular momentum \(l=\)`l`.
      procedure, public :: get_rad_fun
      !> Get the radial functions \(g^{\alpha,\pm}_\tilde{\lambda}(r)\) for a given index \(\tilde{\lambda}=\)`lam`, 
      !> atom index \(\alpha=\)`ias` and angular momentum \(l=\)`l` defining the gradient of basis functions.
      procedure, public :: get_gradient_rad_fun
      !> Get transformation matrix to transform between muffin-tin basis and standard (L)APW+LO basis.
      procedure, public :: get_basis_transform
      !> Transform the eigenvectors \(C^{n{\bf k}}_\mu\) in the standard (L)APW+LO basis 
      !> into the muffin-tin basis representation \(C^{n{\bf k}\alpha}_\lambda\)
      procedure, public :: transform_evec
  end type mt_basis_type
  !> Constructor for [[mt_basis_type(type)]] object.
  interface mt_basis_type
    module procedure :: setup_mt_basis
  end interface

  public :: mt_basis_type, generate_non_zero_clebsch_gordan

contains

  !> Setup the muffin-tin basis object.
  function setup_mt_basis( rad_grid, n_rad_grid, apw_rad_fun, lo_rad_fun, lmax_basis, apword, nlorb, lorbl ) result( this )
    use asserts
    !> radial grids on which the radial functions are given for each species
    real(dp), intent(in) :: rad_grid(:,:)
    !> number of radial grid points for each species
    integer, intent(in) :: n_rad_grid(:)
    !> (L)APW radial functions \(u^\alpha_{l,\xi}(r)\)
    real(dp), intent(in) :: apw_rad_fun(:,:,:,:,:)
    !> LO radial functions \(f^L(r)\)
    real(dp), intent(in) :: lo_rad_fun(:,:,:,:)
    !> maximum angular momentum considered for basis functions
    integer, intent(in) :: lmax_basis
    !> maximum APW linearization order for each angular momentum and species
    integer, intent(in) :: apword(0:,:)
    !> number of local orbitals for all species
    integer, intent(in) :: nlorb(:)
    !> angular momentum for each local orbital and species
    integer, intent(in) :: lorbl(:,:)

    !> `mt_basis_type` type object
    type(mt_basis_type) :: this

    integer :: nspecies, natoms, nradmax, apwordmax, lmaxapw, lmaxlo, nlomax, dordmax
    integer :: is, ir, l, m, lm, o, ilo, lam

    ! defer sizes from input arrays
    nspecies = size( rad_grid, dim=2 )
    natoms = size( apw_rad_fun, dim=5 )
    nradmax = maxval( n_rad_grid )
    lmaxapw = size( apw_rad_fun, dim=4 ) - 1
    apwordmax = maxval( apword(0:lmaxapw, 1:nspecies) )
    nlomax = maxval( nlorb(1:nspecies) )
    lmaxlo = maxval( lorbl(1:nlomax, 1:nspecies) )
    dordmax = size( apw_rad_fun, dim=2 ) - 1

    call assert( lmax_basis >= 0, &
      'Maximum l of basis functions must not be negative.' )
    call assert( size( n_rad_grid ) == nspecies, &
      'Radial grids and number of radial grid points given for different number of species.' )
    call assert( size( lo_rad_fun, dim=4 ) == natoms, &
      'Radial APW functions and LO functions given for different number of atoms.' )
    call assert( size( apw_rad_fun, dim=1 ) >= nradmax, &
      'Radial APW functions not given for enough radial points.' )
    call assert( size( lo_rad_fun, dim=1 ) >= nradmax, &
      'Radial LO functions not given for enough radial points.' )
    call assert( size( lo_rad_fun, dim=2 ) == dordmax + 1, &
      'Radial derivative for APW functions and LO functions given up to different order.' )
    call assert( size( apw_rad_fun, dim=3 ) >= apwordmax, &
      'Radial APW functions not given for large enough linearization orders.' )
    call assert( size( lo_rad_fun, dim=3 ) >= nlomax, &
      'Radial LO functions not given for enough LOs.' )

    call this%destroy

    this%lmax_basis = lmax_basis

    ! assign radial grids
    allocate( this%n_rad_grid(nspecies), source=n_rad_grid )
    allocate( this%rad_grid(nradmax, nspecies), source=0.0_dp )
    do is = 1, nspecies
      do ir = 1, this%n_rad_grid(is)
        this%rad_grid(ir, is) = rad_grid(ir, is)
      end do
    end do

    ! assign radial functions
    allocate( this%apw_rad_fun(nradmax, 0:dordmax, apwordmax, 0:lmaxapw, natoms) )
    allocate( this%lo_rad_fun(nradmax, 0:dordmax, nlomax, natoms) )
    this%apw_rad_fun = apw_rad_fun(1:nradmax, :, 1:apwordmax, :, :)
    this%lo_rad_fun = lo_rad_fun(1:nradmax, :, 1:nlomax, :)

    ! build index maps
    allocate( this%n_rad_fun(0:lmax_basis, nspecies), source=0 )
    allocate( this%rad_fun_to_apw_lo(apwordmax+nlomax, 0:lmax_basis, nspecies), source=0 )
    do is = 1, nspecies
      do l = 0, min( lmaxapw, lmax_basis )
        do o = 1, apword(l, is)
          this%n_rad_fun(l, is) = this%n_rad_fun(l, is) + 1
          this%rad_fun_to_apw_lo(this%n_rad_fun(l, is), l, is) = o
        end do
      end do
      do l = 0, min( lmaxlo, lmax_basis )
        do ilo = 1, nlorb(is)
          if( lorbl(ilo, is) /= l ) cycle
          this%n_rad_fun(l, is) = this%n_rad_fun(l, is) + 1
          this%rad_fun_to_apw_lo(this%n_rad_fun(l, is), l, is) = -ilo
        end do
      end do
    end do
    this%n_rad_fun_max = maxval( this%n_rad_fun )

    allocate( this%n_basis_fun(nspecies), source=0 )
    allocate( this%idx_basis_fun((lmax_basis+1)**2, this%n_rad_fun_max,  nspecies), source=0 )
    lm = 0
    do l = 0, lmax_basis
      do m = -l, l
        lm = lm + 1
        do is = 1, nspecies
          do lam = 1, this%n_rad_fun(l, is)
            this%n_basis_fun(is) = this%n_basis_fun(is) + 1
            this%idx_basis_fun(lm, lam, is) = this%n_basis_fun(is)
          end do
        end do
      end do
    end do
    this%n_basis_fun_max = maxval( this%n_basis_fun )
  end function setup_mt_basis

  subroutine destroy( this )
    !> `mt_basis_type` type object
    class(mt_basis_type), intent(inout) :: this
    this%n_rad_fun_max = 0
    this%n_basis_fun_max = 0
    this%lmax_basis = 0
    if( allocated( this%apw_rad_fun ) ) deallocate( this%apw_rad_fun )
    if( allocated( this%lo_rad_fun ) ) deallocate( this%lo_rad_fun )
    if( allocated( this%rad_grid ) ) deallocate( this%rad_grid )
    if( allocated( this%n_rad_grid ) ) deallocate( this%n_rad_grid )
    if( allocated( this%n_rad_fun ) ) deallocate( this%n_rad_fun )
    if( allocated( this%n_basis_fun ) ) deallocate( this%n_basis_fun )
    if( allocated( this%idx_basis_fun ) ) deallocate( this%idx_basis_fun )
    if( allocated( this%rad_fun_to_apw_lo ) ) deallocate( this%rad_fun_to_apw_lo )
  end subroutine destroy

  !> Get the radial function \(g^\alpha_\tilde{\lambda}(r)\) corresponding to basis functions 
  !> \[ \phi^\alpha_\lambda({\bf r}) = g^\alpha_\tilde{\lambda}(r_\alpha)\, Y_{l_\lambda m_\lambda}(\hat{\bf r}_\alpha) \]
  !> for a given index \(\tilde{\lambda}=\)`lam`, atom index \(\alpha=\)`ias` and angular momentum \(l=\)`l`.
  function get_rad_fun( this, l, is, ias, lam, radial_derivative ) result( rad_fun )
    use asserts
    !> `mt_basis_type` type object
    class(mt_basis_type), intent(in) :: this
    !> angular momentum \(l\)
    integer, intent(in) :: l
    !> species index
    integer, intent(in) :: is
    !> index of atom \(\alpha\)
    integer, intent(in) :: ias
    !> index \(\tilde{\lambda}\) of radial function
    integer, intent(in) :: lam
    !> order of radial derivative (default: `0`)
    integer, optional, intent(in) :: radial_derivative

    !> radial function
    real(dp), allocatable :: rad_fun(:)

    real(dp), allocatable :: f(:), cf(:,:)

    integer :: dord, dordmax, i

    dord = 0
    if( present( radial_derivative ) ) dord = radial_derivative

    call assert( dord >= 0, &
      'Radial derivative order must not be negative.' )
    call assert( lam > 0 .and. lam <= this%n_rad_fun(l, is), &
      'Index `lam` out of bounds.' )

    if( allocated( rad_fun ) ) deallocate( rad_fun )
    allocate( rad_fun(this%n_rad_grid(is)) )

    dordmax = ubound( this%apw_rad_fun, dim=2 )
    i = this%rad_fun_to_apw_lo(lam, l, is)
    call assert( i /= 0, &
      'Unable to assign given combination of `l`, `is` and `lam` to either an APW or an LO radial function' )
    if( i > 0 ) then
      rad_fun = this%apw_rad_fun(1:this%n_rad_grid(is), min(dord, dordmax), i, l, ias)
    else
      rad_fun = this%lo_rad_fun(1:this%n_rad_grid(is), min(dord, dordmax), -i, ias)
    end if

    ! requested derivative order is not available and needs to be computed
    if( dord > dordmax ) then
      allocate( f(this%n_rad_grid(is)), cf(3, this%n_rad_grid(is)) )
      f = rad_fun
      call fderiv( dord-dordmax, this%n_rad_grid(is), this%rad_grid(1:this%n_rad_grid(is), is), f, rad_fun, cf )
      deallocate( f, cf )
    end if
  end function get_rad_fun

  !> The gradient of basis functions can be expressed as
  !> \[
  !> {\bf \nabla}\phi^\alpha_\lambda({\bf r}) = {\bf \nabla}[g^\alpha_{\tilde{\lambda}}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha)] = 
  !> \left[ \frac{l}{2l+1} \right]^\frac{1}{2} \left[ \frac{l+1}{r_\alpha} + \frac{{\rm d}}{{\rm d}r_\alpha} \right]
  !> g^\alpha_{\tilde{\lambda}}(r_\alpha) \, {\bf Y}_{l\,l-1\,m}(\hat{\bf r}_\alpha) +
  !> \left[ \frac{l+1}{2l+1} \right]^\frac{1}{2} \left[ \frac{l}{r_\alpha} - \frac{{\rm d}}{{\rm d}r_\alpha} \right]
  !> g^\alpha_{\tilde{\lambda}}(r_\alpha) \, {\bf Y}_{l\,l+1\,m}(\hat{\bf r}_\alpha) \;,
  !> \]
  !> where \({\bf Y}_{ll'm}(\hat{\bf r}_\alpha)\) are the so called vector spherical harmonics
  !> \[
  !> {\bf Y}_{ll'm}(\hat{\bf r}_\alpha) = \sum_{m'=-l'}^{l'} Y_{l'm'}(\hat{\bf r}_\alpha) \begin{pmatrix}
  !> [C(l'\,1\,l|m'\,-1\,m) - C(l'\,1\,l|m'\,1\,m)]/\sqrt{2} \\
  !> -{\rm i}[C(l'\,1\,l|m'\,-1\,m) + C(l'\,1\,l|m'\,1\,m)]/\sqrt{2} \\
  !> C(l'\,1\,l|m'\,0\,m) \end{pmatrix} \;,
  !> \]
  !> with \(C(l_1\,l_2\,l_3|m_1\,m_2\,m_3)\) being Clebsch-Gordan coefficients
  !> (see [[wigner3j_symbol(module):clebsch_gordan(function)]]). In a simplified notation, we
  !> find
  !> \[
  !> {\bf \nabla}[g^\alpha_{\tilde{\lambda}}(r) \, Y_{lm}(\hat{\bf r})] = \sum_{\pm}
  !> g^{\alpha,\pm}_{\tilde{\lambda}}(r) \, \sum_{m'=-(l\pm 1)}^{l\pm1} {\bf c}^{\pm}_{mm',l} \, Y_{l\pm1 m'}(\hat{\bf r}) \;,
  !> \]
  !> with 
  !> \[ g^{\alpha,-}_\tilde{\lambda}(r) = \left[ \frac{l}{2l+1} \right]^\frac{1}{2} \left[ \frac{l+1}{r} + \frac{{\rm d}}{{\rm d}r} \right] g^\alpha_\tilde{\lambda}(r) \; , \]
  !> and 
  !> \[ g^{\alpha,+}_\tilde{\lambda}(r) = \left[ \frac{l+1}{2l+1} \right]^\frac{1}{2} \left[ \frac{l}{r} - \frac{{\rm d}}{{\rm d}r} \right] g^\alpha_\tilde{\lambda}(r) \; , \]
  !> and
  !> \[ {\bf c}^{\pm}_{mm',l} = \begin{pmatrix}
  !> [C(l \pm 1\,1\,l|m'\,-1\,m) - C(l \pm 1\,1\,l|m'\,1\,m)]/\sqrt{2} \\
  !> {\rm i}[C(l \pm 1\,1\,l|m'\,-1\,m) + C(l \pm 1\,1\,l|m'\,1\,m)]/\sqrt{2} \\
  !> C(l \pm 1\,1\,l|m'\,0\,m) \end{pmatrix} \;. \]
  !> This routine returns \(g^{\alpha,\pm}_\tilde{\lambda}(r)\) or one of its radial derivatives
  !> depending on the input value or the arguments `grad` and `radial_derivative`.
  function get_gradient_rad_fun( this, l, is, ias, lam, grad, radial_derivative ) result( rad_fun )
    use asserts
    !> `mt_basis_type` type object
    class(mt_basis_type), intent(in) :: this
    !> angular momentum \(l\)
    integer, intent(in) :: l
    !> species index of atom \(\alpha\)
    integer, intent(in) :: is
    !> index of atom \(\alpha\)
    integer, intent(in) :: ias
    !> index \(\tilde{\lambda}\) of radial function
    integer, intent(in) :: lam
    !> gradient term 
    !> (`>0` for \(g^{\alpha,+}_\tilde{\lambda}\), `<0` for \(g^{\alpha,-}_\tilde{\lambda}\), `=0` for \(g^\alpha_\tilde{\lambda}\) )
    integer, intent(in) :: grad
    !> order of radial derivative (default: `0`)
    integer, optional, intent(in) :: radial_derivative

    !> radial function
    real(dp), allocatable :: rad_fun(:)

    real(dp), allocatable :: f(:), cf(:,:)

    integer :: dord, dordmax, i
    real(dp) :: lfac

    dord = 0
    if( present( radial_derivative ) ) dord = radial_derivative

    call assert( dord >= 0, &
      'Radial derivative order must not be negative.' )
    call assert( lam > 0 .and. lam <= this%n_rad_fun(l, is), &
      'Index `lam` out of bounds.' )

    if( grad == 0 ) then
      rad_fun = this%get_rad_fun( l, is, ias, lam, radial_derivative=dord )
      return
    else if( grad < 0 ) then
      lfac = sqrt( dble( l ) / dble( 2 * l + 1 ) )
    else
      lfac = sqrt( dble( l + 1 ) / dble( 2 * l + 1 ) )
    end if

    if( allocated( rad_fun ) ) deallocate( rad_fun )
    allocate( rad_fun(this%n_rad_grid(is)) )

    dordmax = ubound( this%apw_rad_fun, dim=2 )
    i = this%rad_fun_to_apw_lo(lam, l, is)
    call assert( i /= 0, &
      'Unable to assign given combination of `l`, `is` and `lam` to either an APW or an LO radial function' )
    if( i > 0 ) then
      if( grad < 0 ) then
        rad_fun = lfac * ( (l + 1) * this%apw_rad_fun(1:this%n_rad_grid(is), 0, i, l, ias) / this%rad_grid(1:this%n_rad_grid(is), is) &
                          + this%apw_rad_fun(1:this%n_rad_grid(is), 1, i, l, ias) )
      else
        rad_fun = lfac * (       l * this%apw_rad_fun(1:this%n_rad_grid(is), 0, i, l, ias) / this%rad_grid(1:this%n_rad_grid(is), is) &
                          - this%apw_rad_fun(1:this%n_rad_grid(is), 1, i, l, ias) )
      end if
    else
      if( grad < 0 ) then
        rad_fun = lfac * ( (l + 1) * this%lo_rad_fun(1:this%n_rad_grid(is), 0, -i, ias) / this%rad_grid(1:this%n_rad_grid(is), is) &
                          + this%lo_rad_fun(1:this%n_rad_grid(is), 1, -i, ias) )
      else
        rad_fun = lfac * (       l * this%lo_rad_fun(1:this%n_rad_grid(is), 0, -i, ias) / this%rad_grid(1:this%n_rad_grid(is), is) &
                          - this%lo_rad_fun(1:this%n_rad_grid(is), 1, -i, ias) )
      end if
    end if

    ! compute radial derivative if requested
    if( dord > 0 ) then
      allocate( f(this%n_rad_grid(is)), cf(3, this%n_rad_grid(is)) )
      f = rad_fun
      call fderiv( dord, this%n_rad_grid(is), this%rad_grid(1:this%n_rad_grid(is), is), f, rad_fun, cf )
      deallocate( f, cf )
    end if
  end function get_gradient_rad_fun

  !> Get the transformation matrix \(T^\alpha_{\lambda\mu}\) to transform any vector \(w^\alpha_\lambda\) 
  !> given in the muffin-tin basis indices \(\lambda\) into a vector \(v^\alpha_\mu\) given in the
  !> standard (L)APW+LO basis indices
  !> \[ \mu = \begin{cases}{\bf G+k} & \text{ for (L)APWs} \\ L & \text{ for LOs} \end{cases} \]
  !> such that
  !> \[ {\bf v}^\alpha = {{\bf T}^\alpha}^\top \cdot {\bf w}^\alpha \; . \]
  subroutine get_basis_transform( this, is, ngk, nlotot, idxlo, apwalm, transform, use_local_orbitals )
    use constants, only: zzero, zone
    !> `mt_basis_type` type object
    class(mt_basis_type), intent(in) :: this
    !> species index of atom \(\alpha\)
    integer, intent(in) :: is
    !> number of \({\bf G+k}\) vectors
    integer, intent(in) :: ngk
    !> total number of local orbitals
    integer, intent(in) :: nlotot
    !> index of the LO \(f^\alpha_L(r)\, Y_{lm}(\hat{\bf r})\) in the standard basis \(\mu\)
    !> for all \((l,m)=\)`lm` and \(L=\)`ilo`
    integer, intent(in) :: idxlo(:,:)
    !> (L)APW matching coefficients \(A^\alpha_{lm,\xi,{\bf G+k}}\)
    complex(dp), intent(in) :: apwalm(:,:,:)
    !> transformation matrix \({\bf T}^\alpha\)
    complex(dp), allocatable, intent(out) :: transform(:,:)
    !> include local orbitals or not (default: `.true.`)
    logical, optional, intent(in) :: use_local_orbitals

    integer :: l, lm, lam, i
    logical :: uselo

    uselo = .true.
    if( present( use_local_orbitals ) ) uselo = use_local_orbitals

    if( allocated( transform ) ) deallocate( transform )
    allocate( transform(this%n_basis_fun(is), ngk+nlotot), source=zzero )

!$omp parallel default( shared ) private( l, lm, lam, i )
!$omp do
    do lm = 1, (this%lmax_basis+1)**2
      l = ceiling( sqrt( dble( lm ) ) - (1.0_dp + 1e-16_dp) )
      do lam = 1, this%n_rad_fun(l, is)
        i = this%rad_fun_to_apw_lo(lam, l, is)
        if( i > 0 ) then
          transform(this%idx_basis_fun(lm, lam, is), 1:ngk) = apwalm(1:ngk, i, lm)
        else if( i < 0 .and. uselo ) then
          transform(this%idx_basis_fun(lm, lam, is), ngk+idxlo(lm, -i)) = zone
        end if
      end do
    end do
!$omp end do
!$omp end parallel
  end subroutine get_basis_transform

  !> Transform the eigenvectors \(C^{n{\bf k}}_\mu\) given in the standard (L)APW+LO basis indices
  !> \[ \mu = \begin{cases}{\bf G+k} & \text{ for (L)APWs} \\ L & \text{ for LOs} \end{cases} \]
  !> into the muffin-tin basis representation 
  !> \[ C^{n{\bf k}\alpha}_\lambda = \begin{cases} 
  !>    \sum_{\bf G+k} C^{n{\bf k}}_{\bf G+k}\, A^\alpha_{lm, \xi, {\bf G+k}} & \lambda \rightarrow \text{ (L)APW} \\
  !>    C^{n{\bf k}}_L\, \delta_{\alpha \alpha_L} & \lambda \rightarrow \text{ LO} \; . \end{cases} \]
  subroutine transform_evec( this, is, ngk, nlotot, idxlo, apwalm, evec, evec_mt, use_local_orbitals )
    use asserts
    use constants, only: zzero, zone
    !> `mt_basis_type` type object
    class(mt_basis_type), intent(in) :: this
    !> species index of atom \(\alpha\)
    integer, intent(in) :: is
    !> number of \({\bf G+k}\) vectors
    integer, intent(in) :: ngk
    !> total number of local orbitals
    integer, intent(in) :: nlotot
    !> index of the LO \(f^\alpha_L(r)\, Y_{lm}(\hat{\bf r})\) in the standard basis \(\mu\)
    !> for all \((l,m)=\)`lm` and \(L=\)`ilo`
    integer, intent(in) :: idxlo(:,:)
    !> (L)APW matching coefficients \(A^\alpha_{lm,\xi,{\bf G+k}}\)
    complex(dp), intent(in) :: apwalm(:,:,:)
    !> eigenvectors \(C^{n{\bf k}}_\mu\) in standard basis
    complex(dp), intent(in) :: evec(:,:)
    !> eigenvectors \(C^{n{\bf k}\alpha}_\lambda\) in muffin-tin basis
    complex(dp), allocatable, intent(out) :: evec_mt(:,:)
    !> include local orbitals or not (default: `.true.`)
    logical, optional, intent(in) :: use_local_orbitals

    integer :: l, lm, lam, i, adim(3), vdim(2)
    logical :: uselo

    complex(dp), allocatable :: contracted(:,:,:)

    uselo = .true.
    if( present( use_local_orbitals ) ) uselo = use_local_orbitals

    adim = shape( apwalm )
    vdim = shape( evec )

    call assert( vdim(1) >= ngk + nlotot, &
      'First dimension of eigenvector array too small.' )

    if( allocated( evec_mt ) ) deallocate( evec_mt )
    allocate( evec_mt(this%n_basis_fun(is), vdim(2)), source=zzero )

    allocate( contracted(vdim(2), adim(2), adim(3)) )
    call zgemm( 't', 'n', vdim(2), adim(2)*adim(3), ngk, zone, &
      evec, vdim(1), &
      apwalm, adim(1), zzero, &
      contracted, vdim(2) )

!$omp parallel default( shared ) private( l, lm, lam, i )
!$omp do
    do lm = 1, (this%lmax_basis+1)**2
      l = ceiling( sqrt( dble( lm ) ) - (1.0_dp + 1e-16_dp) )
      do lam = 1, this%n_rad_fun(l, is)
        i = this%rad_fun_to_apw_lo(lam, l, is)
        if( i > 0 ) then
          evec_mt(this%idx_basis_fun(lm, lam, is), :) = contracted(:, i, lm)
        else if( i < 0 .and. uselo ) then
          evec_mt(this%idx_basis_fun(lm, lam, is), :) = evec(ngk+idxlo(lm, -i), :)
        end if
      end do
    end do
!$omp end do
!$omp end parallel

    deallocate( contracted )
  end subroutine transform_evec

  subroutine generate_non_zero_clebsch_gordan( lmax, eps, cg_num, cg_lm, cg_val )
    use constants, only: zzero, zone, sqrt_two
    use wigner3j_symbol, only: clebsch_gordan
    integer, intent(in) :: lmax
    real(dp), intent(in) :: eps
    integer, allocatable, intent(out) :: cg_num(:,:,:)
    integer, allocatable, intent(out) :: cg_lm(:,:,:,:)
    complex(dp), allocatable, intent(out) :: cg_val(:,:,:,:)

    integer :: i, l, m, lm, ll, mm, lmmax, sig
    real(dp) :: cg(-1:1), t1

    lmmax = (lmax + 1)**2

    if( allocated( cg_num ) ) deallocate( cg_num )
    if( allocated( cg_lm ) ) deallocate( cg_lm )
    if( allocated( cg_val ) ) deallocate( cg_val )
    allocate( cg_num(0:3, -1:1, lmmax), source=0 )
    allocate( cg_lm((lmax + 2)**2, 0:3, -1:1, lmmax), source=0 )
    allocate( cg_val((lmax + 2)**2, 0:3, -1:1, lmmax), source=zzero )

    do sig = -1, 1
      do l = max(0, -sig), lmax
        do m = -l, l
          lm = l * (l + 1) + m + 1
          if( sig == 0 ) then
            cg_num(:, 0, lm) = 1
            cg_lm(1, :, 0, lm) = lm
            cg_val(1, :, 0, lm) = zone
          else
            ll = l + sig
            do mm = -ll, ll
              cg = [(clebsch_gordan( ll, 1, l, mm, i, m ), i=-1, 1)]
              ! x-component
              t1 = cg(-1) - cg(1)
              if( abs( t1 ) > eps ) then
                cg_num(1, sig, lm) = cg_num(1, sig, lm) + 1
                cg_lm(cg_num(1, sig, lm), 1, sig, lm) = ll * (ll + 1) + mm + 1
                cg_val(cg_num(1, sig, lm), 1, sig, lm) = cmplx( t1 / sqrt_two, 0, dp )
              end if
              ! y-component
              t1 = cg(-1) + cg(1)
              if( abs( t1 ) > eps ) then
                cg_num(2, sig, lm) = cg_num(2, sig, lm) + 1
                cg_lm(cg_num(2, sig, lm), 2, sig, lm) = ll * (ll + 1) + mm + 1
                cg_val(cg_num(2, sig, lm), 2, sig, lm) = cmplx( 0, - t1 / sqrt_two, dp )
              end if
              ! z-component
              t1 = cg(0)
              if( abs( t1 ) > eps ) then
                cg_num(3, sig, lm) = cg_num(3, sig, lm) + 1
                cg_lm(cg_num(3, sig, lm), 3, sig, lm) = ll * (ll + 1) + mm + 1
                cg_val(cg_num(3, sig, lm), 3, sig, lm) = cmplx( t1, 0, dp )
              end if
            end do
          end if
        end do
      end do
    end do
  end subroutine generate_non_zero_clebsch_gordan

end module muffin_tin_basis
