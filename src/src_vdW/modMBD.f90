!> Many-body dispersion corrections under periodic boundary conditions.
!>
!> This module implements the range-separated self-consistent screening
!> variant of the many-body dispersion method (MBD@rsSCS). It represents each
!> atom by an isotropic quantum harmonic oscillator and treats their coupled
!> dipole response. The method separates the short-range screening already
!> described by the density functional from the long-range correlation energy.
!>
!> The implementation proceeds in three stages. First, the Hirshfeld volume
!> ratio \(\nu_i=V_i/V_i^{\mathrm{free}}\) scales the free-atom TS parameters,
!> \[
!>   \alpha_i^{\mathrm{TS}}=\nu_i\alpha_i^{\mathrm{free}}, \qquad
!>   R_{0,i}^{\mathrm{TS}}=\nu_i^{1/3}R_{0,i}^{\mathrm{free}}.
!> \]
!> The dynamic TS polarizability is
!> \[
!>   \alpha_i^{\mathrm{TS}}(\mathrm{i}\omega)=
!>   \frac{\alpha_i^{\mathrm{TS}}}
!>        {1+(\omega/\omega_i)^2}, \qquad
!>   \omega_i=\frac{4C_{6,i}}{3(\alpha_i^{\mathrm{free}})^2}.
!> \]
!> Second, `shortMBD` obtains the screened response from the inverse of the
!> short-range matrix
!> \[
!>   B^{\mathrm{SR}}_{i\alpha,j\beta}=
!>   \frac{\delta_{ij}\delta_{\alpha\beta}}
!>        {\alpha_i^{\mathrm{TS}}(\mathrm{i}\omega)}+
!>   T^{\mathrm{SR}}_{i\alpha,j\beta}, \qquad
!>   \bar\alpha_i=\frac{1}{3}\operatorname{Tr}
!>   \sum_j\left(B^{\mathrm{SR}}\right)^{-1}_{ij}.
!> \]
!> Here \(T^{\mathrm{SR}}\) is the Gaussian-attenuated, short-range complement
!> of the dipole interaction. The bare tensor used in the lattice sums is
!> \(T_{\alpha\beta}(\mathbf{R})=(3R_\alpha R_\beta-
!> R^2\delta_{\alpha\beta})/R^5\).
!>
!> Finally, `longMBD` constructs the range-separated long-range tensor at each
!> Brillouin-zone point,
!> \[
!>   T^{\mathrm{LR}}_{ij}(\mathbf{k})=
!>   \sum_{\mathbf{L}} f_{\mathrm{damp}}(R_{ij\mathbf{L}})
!>   T(\mathbf{R}_{ij\mathbf{L}})
!>   e^{-\mathrm{i}\mathbf{k}\cdot\mathbf{L}},
!> \]
!> and evaluates the imaginary-frequency correlation energy
!> \[
!>   E_{\mathrm{MBD}}=\frac{1}{2\pi}\sum_{\mathbf{k}}w_{\mathbf{k}}
!>   \int_0^\infty\!d\omega\left[
!>   \ln\det\left(\mathbf{I}-\mathbf{A}T^{\mathrm{LR}}(\mathbf{k})\right)
!>   +\operatorname{Tr}\left(\mathbf{A}T^{\mathrm{LR}}(\mathbf{k})\right)
!>   \right],
!> \]
!> where \(\mathbf{A}\) is diagonal in atomic and Cartesian indices with entries
!> \(\bar\alpha_i(\mathrm{i}\omega)\). `longfMBD` evaluates
!> \(\mathbf{F}_I=-\partial E_{\mathrm{MBD}}/\partial\mathbf{R}_I\), including
!> derivatives of both the screened polarizabilities and the range-separated
!> interaction tensor.
!>
!> ### References
!>
!> 1. A. Tkatchenko, R. A. DiStasio Jr., R. Car, and M. Scheffler,
!>    [Accurate and Efficient Method for Many-Body van der Waals Interactions,
!>    *Phys. Rev. Lett.* **108**, 236402 (2012)]
!>    (https://doi.org/10.1103/PhysRevLett.108.236402).
!> 2. A. Ambrosetti, A. M. Reilly, R. A. DiStasio Jr., and A. Tkatchenko,
!>    [Long-range correlation energy calculated from coupled atomic response
!>    functions, *J. Chem. Phys.* **140**, 18A508 (2014)]
!>    (https://doi.org/10.1063/1.4865104).
!> 3. T. Bučko, S. Lebègue, T. Gould, and J. G. Ángyán,
!>    [Many-body dispersion corrections for periodic systems: an efficient
!>    reciprocal space implementation, *J. Phys.: Condens. Matter* **28**,
!>    045201 (2016)] (https://doi.org/10.1088/0953-8984/28/4/045201).
module modMBD
  use precision, only: dp

  Implicit none
  private
  public :: longMBD, longfMBD

  integer, parameter :: omegadim=19
  real(dp), allocatable :: omg(:)
  real(dp), allocatable :: omegaweight(:)
  real(dp), allocatable :: xyz(:,:)
  !> maximal interaction distance
  real(dp) :: rsrcutoff

 contains

  !> Determine lattice replication counts from the MBD cutoff radius.
  Subroutine getlatttrans(lattrans)
    Use modinput, only: input
    Use vdw_general_routines, only: cross
    Implicit None
    !> Number of translated unit cells in each lattice direction.
    Integer, Intent(out) :: lattrans(3)
    real(dp) :: vec(3,3)
    Integer :: icount

    vec(:,1) = cross(input%structure%crystal%basevect(:,2),input%structure%crystal%basevect(:,3))
    vec(:,2) = cross(input%structure%crystal%basevect(:,1),input%structure%crystal%basevect(:,3))
    vec(:,3) = cross(input%structure%crystal%basevect(:,1),input%structure%crystal%basevect(:,2))
    Do icount = 1,3
       vec(:,icount) = vec(:,icount)/Sqrt(Dot_product(vec(:,icount),vec(:,icount)))!normalize
       lattrans(icount) = Int(abs(rsrcutoff/(Dot_product(input%structure%crystal%basevect(:,icount),vec(:,icount)))))+1
    End Do
  End Subroutine

  !> Diagonal Cartesian component of the short-range interaction-tensor derivative.
  !>
  !> Evaluates the Gaussian-attenuated derivative for three equal Cartesian
  !> indices, \(\partial T_{aa}/\partial r_a\).
  real(dp)  function derivTsr_a(a,rijn,sigij)
     use modmain, only: pi

     implicit none
     !> Cartesian component of the interatomic vector.
     real(dp), intent(in) :: a
     !> Distance between atoms \(i\) and \(j\).
     real(dp), intent(in) :: rijn
     !> Gaussian attenuation length.
     real(dp), intent(in) :: sigij

     real(dp) :: t1, t2, t3
     real(dp) :: rsij
     real(dp) :: rijn2, rijn4
     real(dp) :: rijn_inv, rijn_inv2, rijn_inv3, rijn_inv4, rijn_inv5, rijn_inv6, rijn_inv7

     rsij=rijn/sigij
     rijn2=rijn*rijn
     rijn4=rijn2*rijn2
     rijn_inv=1.0_dp/rijn
     rijn_inv2=rijn_inv*rijn_inv
     rijn_inv3=rijn_inv2*rijn_inv
     rijn_inv4=rijn_inv2*rijn_inv2
     rijn_inv5=rijn_inv4*rijn_inv
     rijn_inv6=rijn_inv3*rijn_inv3
     rijn_inv7=rijn_inv6*rijn_inv

     t1=(-9.0_dp*a*rijn_inv5+15.0_dp*(a**3)*rijn_inv7)*erf(rsij)
     t2=3.0_dp/sqrt(pi)*a*rijn_inv4*1.0_dp/(sigij**3)*(6.0_dp*(sigij**2)+4*rijn2)*exp(-(rsij**2))
     t3=-1.0_dp/sqrt(pi)*(a**3)*rijn_inv6*1.0_dp/(sigij**5)*(30.0_dp*(sigij**4)+20.0_dp*(sigij**2)*rijn2+8.0_dp*rijn4)*exp(-(rsij**2))
     derivTsr_a=(t1+t2+t3)

  end function

  !> Two-equal-components contribution to the short-range interaction-tensor derivative.
  !>
  !> Evaluates \(\partial T_{bb}/\partial r_a\) for \(a \ne b\), including
  !> the error-function attenuation and its Gaussian derivative.
  real(dp)  function derivTsr_ab(a,b,rijn,sigij)
     use modmain, only: pi

     implicit none
     !> Single Cartesian component of the interatomic vector.
     real(dp), intent(in) :: a
     !> Repeated Cartesian component of the interatomic vector.
     real(dp), intent(in) :: b
     !> Distance between atoms \(i\) and \(j\).
     real(dp), intent(in) :: rijn
     !> Gaussian attenuation length.
     real(dp), intent(in) :: sigij

     real(dp) :: t1, t2, t3
     real(dp) :: rsij
     real(dp) :: rijn2, rijn4
     real(dp) :: rijn_inv, rijn_inv2, rijn_inv3, rijn_inv4, rijn_inv5, rijn_inv6, rijn_inv7

     rsij=rijn/sigij
     rijn2=rijn*rijn
     rijn4=rijn2*rijn2
     rijn_inv=1.0_dp/rijn
     rijn_inv2=rijn_inv*rijn_inv
     rijn_inv3=rijn_inv2*rijn_inv
     rijn_inv4=rijn_inv2*rijn_inv2
     rijn_inv5=rijn_inv4*rijn_inv
     rijn_inv6=rijn_inv3*rijn_inv3
     rijn_inv7=rijn_inv6*rijn_inv

     t1=-(3.0_dp*a*rijn_inv5-15.0_dp*a*(b**2)*rijn_inv7)*erf(rsij)
     t2=1.0_dp/sqrt(pi)*a*rijn_inv4*1.0_dp/(sigij**3)*(6.0_dp*(sigij**2)+4*rijn2)*exp(-(rsij**2))
     t3=-1.0_dp/sqrt(pi)*a*(b**2)*rijn_inv6*1.0_dp/(sigij**5)*(30.0_dp*(sigij**4)+20.0_dp*(sigij**2)*rijn2+8.0_dp*rijn4)*exp(-(rsij**2))
     derivTsr_ab=(t1+t2+t3)

  end function

  !> Three-distinct-components contribution to the short-range interaction-tensor derivative.
  !>
  !> Evaluates \(\partial T_{bc}/\partial r_a\) for three distinct Cartesian
  !> indices using the Gaussian-attenuated dipole interaction.
  real(dp)  function derivTsr_abc(a,b,c,rijn,sigij)
     use modmain, only: pi

     implicit none
     !> First Cartesian component of the interatomic vector.
     real(dp), intent(in) :: a
     !> Second Cartesian component of the interatomic vector.
     real(dp), intent(in) :: b
     !> Third Cartesian component of the interatomic vector.
     real(dp), intent(in) :: c
     !> Distance between atoms \(i\) and \(j\).
     real(dp), intent(in) :: rijn
     !> Gaussian attenuation length.
     real(dp), intent(in) :: sigij

     real(dp) :: t1, t2, t3
     real(dp) :: rsij
     real(dp) :: rijn2, rijn4
     real(dp) :: rijn_inv, rijn_inv2, rijn_inv3, rijn_inv4, rijn_inv6, rijn_inv7

     rsij=rijn/sigij
     rijn2=rijn*rijn
     rijn4=rijn2*rijn2
     rijn_inv=1.0_dp/rijn
     rijn_inv2=rijn_inv*rijn_inv
     rijn_inv3=rijn_inv2*rijn_inv
     rijn_inv4=rijn_inv2*rijn_inv2
     rijn_inv6=rijn_inv3*rijn_inv3
     rijn_inv7=rijn_inv6*rijn_inv

     t1=(15.0_dp*a*b*c*rijn_inv7)*erf(rsij)
     t2=-1.0_dp/sqrt(pi)*a*b*c*rijn_inv6*1.0_dp/(sigij**5)*(30.0_dp*(sigij**4)+20.0_dp*(sigij**2)*rijn2+8.0_dp*rijn4)*exp(-(rsij**2))
     derivTsr_abc=(t1+t2)

  end function

  !> Diagonal Cartesian component of the long-range interaction-tensor derivative.
  !>
  !> Evaluates \(\partial T_{aa}/\partial r_a\) for the bare dipole tensor
  !> \(T_{ab}=(3r_a r_b-r^2\delta_{ab})/r^5\).
  real(dp)  function derivTlr_a(a,rijn)

     implicit none
     !> Cartesian component of the interatomic vector.
     real(dp), intent(in) :: a
     !> Distance between atoms \(i\) and \(j\).
     real(dp), intent(in) :: rijn

     real(dp) :: t1
     real(dp) :: rijn2, rijn_inv, rijn_inv2, rijn_inv4, rijn_inv7

     rijn2=rijn*rijn
     rijn_inv=1.0_dp/rijn
     rijn_inv2=rijn_inv*rijn_inv
     rijn_inv4=rijn_inv2*rijn_inv2
     rijn_inv7=rijn_inv4*rijn_inv2*rijn_inv

     t1=(-15.0_dp*(a**3)+3.0_dp*rijn2*(3*a))*rijn_inv7
     derivTlr_a=t1

  end function

  !> Two-equal-components contribution to the long-range interaction-tensor derivative.
  !>
  !> Evaluates \(\partial T_{bb}/\partial r_a\) for \(a \ne b\) from the
  !> bare dipole tensor.
  real(dp)  function derivTlr_ab(a,b,rijn)

     implicit none
     !> Single Cartesian component of the interatomic vector.
     real(dp), intent(in) :: a
     !> Repeated Cartesian component of the interatomic vector.
     real(dp), intent(in) :: b
     !> Distance between atoms \(i\) and \(j\).
     real(dp), intent(in) :: rijn

     real(dp) :: t1
     real(dp) :: rijn2, rijn_inv, rijn_inv2, rijn_inv4, rijn_inv7

     rijn2=rijn*rijn
     rijn_inv=1.0_dp/rijn
     rijn_inv2=rijn_inv*rijn_inv
     rijn_inv4=rijn_inv2*rijn_inv2
     rijn_inv7=rijn_inv4*rijn_inv2*rijn_inv

     t1=(-15.0_dp*a*(b**2)+3.0_dp*rijn2*a)*rijn_inv7
     derivTlr_ab=t1

  end function

  !> Three-distinct-components contribution to the long-range interaction-tensor derivative.
  !>
  !> Evaluates \(\partial T_{bc}/\partial r_a\) for three distinct Cartesian
  !> indices from the bare dipole tensor.
  real(dp)  function derivTlr_abc(a,b,c,rijn)

     implicit none
     !> First Cartesian component of the interatomic vector.
     real(dp), intent(in) :: a
     !> Second Cartesian component of the interatomic vector.
     real(dp), intent(in) :: b
     !> Third Cartesian component of the interatomic vector.
     real(dp), intent(in) :: c
     !> Distance between atoms \(i\) and \(j\).
     real(dp), intent(in) :: rijn

     real(dp) :: t1
     real(dp) :: rijn_inv, rijn_inv2, rijn_inv4, rijn_inv7

     rijn_inv=1.0_dp/rijn
     rijn_inv2=rijn_inv*rijn_inv
     rijn_inv4=rijn_inv2*rijn_inv2
     rijn_inv7=rijn_inv4*rijn_inv2*rijn_inv

     t1=-15.0_dp*a*b*c*rijn_inv7
     derivTlr_abc=t1

  end function

  !> Generate the MBD Gauss-Legendre quadrature frequency grid.
  subroutine frequency_grid(n, omega, omega_weight)

     implicit none
     !> Number of quadrature points.
     integer, intent(in) :: n
     !> Quadrature frequencies.
     real(dp), intent(out) :: omega(n)
     !> Quadrature weights.
     real(dp), intent(out) :: omega_weight(n)

     integer :: i
     real(dp), allocatable :: u(:), wu(:)

     omega = 0.0_dp
     omega_weight = 0.0_dp

     allocate(u(n), wu(n))
     call gauleg(0.0_dp, 1.0_dp, u, wu, n)

     do i = n, 1, -1
        omega(i) = u(n-i+1) / (1.0_dp - u(n-i+1))
        omega_weight(i) = wu(n-i+1) * (1.0_dp - u(n-i+1))**(-2)
     end do

     deallocate(u, wu)

  end subroutine

  !> Calculate the long-range interaction tensor.
  !>
  !> Sums the damped dipole tensor over lattice translations inside the real-space
  !> cutoff and applies the Bloch phase \(e^{-i\mathbf{k}\cdot\mathbf{R}}\).
  subroutine tlrij(ip,jp,kvkl,lattrans,damping_const,Svdw,TLR)
    use constants, only: twopi
    use modinput, only: input

    implicit None

  ! arguments
        !> Fractional coordinates of atom \(i\).
        real(dp), intent(in) :: ip(3)
        !> Fractional coordinates of atom \(j\).
        real(dp), intent(in) :: jp(3)
        !> Fractional k-point vector.
        real(dp), intent(in) :: kvkl(3)
        !> Number of translated unit cells in each lattice direction.
        integer, intent(in) :: lattrans(3)
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> Range-separation radius from screened van der Waals radii.
        real(dp), intent(in) :: Svdw
        !> Long-range interaction tensor.
        complex(dp), intent(out) :: TLR(3,3)
  ! local variables
        integer :: t1, t2, t3, al, bt
        real(dp) :: kx, ky, kz, ltransn, rijn, kl, zt1
        real(dp) :: rijn2, rijn_inv, rijn_inv2, rijn_inv4, rijn_inv5
        integer :: shft(3)
        real(dp) :: ij(3), rL0(3), ltrans(3), rij(3)
        real(dp) :: damp
        complex(dp) :: T(3,3)


    kx=twopi*kvkl(1)
    ky=twopi*kvkl(2)
    kz=twopi*kvkl(3)

    ij=ip-jp
    shft = -nint(ij)
    ij = ij + shft

    rL0=ij(1)*input%structure%crystal%basevect(:,1)+ij(2)*input%structure%crystal%basevect(:,2)&
        &+ij(3)*input%structure%crystal%basevect(:,3)

    TLR=0.0_dp
    T=0.0_dp
    do t1=-lattrans(1),lattrans(1)
       do t2=-lattrans(2),lattrans(2)
          do t3=-lattrans(3),lattrans(3)
             ltrans=t1*(input%structure%crystal%basevect(:,1))+t2*(input%structure%crystal%basevect(:,2))&
             &+t3*(input%structure%crystal%basevect(:,3))
             ltransn=sqrt(dot_product(ltrans,ltrans))
             if (ltransn<=rsrcutoff) then
                rij=rL0+ltrans
                rijn=sqrt(dot_product(rij,rij))
                if (rijn>0.0_dp) then
                   rijn2=rijn*rijn
                   rijn_inv=1.0_dp/rijn
                   rijn_inv2=rijn_inv*rijn_inv
                   rijn_inv4=rijn_inv2*rijn_inv2
                   rijn_inv5=rijn_inv4*rijn_inv
                   damp=1.0_dp/(1.0_dp+exp(-damping_const*((rijn/Svdw)-1.0_dp)))
                   kl=kx*(t1+shft(1))+ky*(t2+shft(2))+kz*(t3+shft(3))
                   T=0.0_dp
                   do al=1,3
                      do bt=1,3
                         if (al==bt) then
                            zt1=(3.0_dp*rij(al)*rij(al)-rijn2)*rijn_inv5
                         else
                            zt1=3.0_dp*rij(al)*rij(bt)*rijn_inv5
                         end if
                         T(al,bt)=zt1
                      enddo
                   enddo
                   TLR=TLR+T*damp*exp(cmplx(0.0_dp, -kl, kind=dp))
                end if
             end if
          enddo
       enddo
    enddo

  end subroutine

  !> Calculate the short-range interaction tensor.
  !>
  !> Sums the Gaussian-attenuated dipole tensor over lattice translations and
  !> retains the short-range complement \(1-f_{\mathrm{damp}}\).
  subroutine tsrij(ip,jp,lattrans,sigij,damping_const,Svdw0,TSR)
    use modmain, only: pi
    use modinput, only: input

    implicit None

  ! arguments
        !> Fractional coordinates of atom \(i\).
        real(dp), intent(in) :: ip(3)
        !> Fractional coordinates of atom \(j\).
        real(dp), intent(in) :: jp(3)
        !> Number of translated unit cells in each lattice direction.
        integer, intent(in) :: lattrans(3)
        !> Gaussian attenuation length.
        real(dp), intent(in) :: sigij
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> Range-separation radius from unscreened atoms-in-molecule van der Waals radii.
        real(dp), intent(in) :: Svdw0
        !> Short-range interaction tensor.
        real(dp), intent(out) :: TSR(3,3)
  ! local variables
        integer :: t1, t2, t3, al, bt
        real(dp) :: ij(3), rL0(3), ltrans(3), rij(3), T(3,3)
        real(dp) :: ltransn, rijn, rsij, zt1, zt2, zt3
        real(dp) :: rijn2, rijn_inv, rijn_inv2, rijn_inv4, rijn_inv5
        real(dp) :: damp


    ij=ip-jp
    ij = ij - nint(ij)

    rL0=ij(1)*input%structure%crystal%basevect(:,1)+ij(2)*input%structure%crystal%basevect(:,2)&
        &+ij(3)*input%structure%crystal%basevect(:,3)

    TSR=0.0_dp
    T=0.0_dp
    do t1=-lattrans(1),lattrans(1)
       do t2=-lattrans(2),lattrans(2)
          do t3=-lattrans(3),lattrans(3)
             ltrans=t1*(input%structure%crystal%basevect(:,1))+t2*(input%structure%crystal%basevect(:,2))&
             &+t3*(input%structure%crystal%basevect(:,3))
             ltransn=sqrt(dot_product(ltrans,ltrans))
             if (ltransn<=rsrcutoff) then
                rij=rL0+ltrans
                rijn=sqrt(dot_product(rij,rij))
                if (rijn>0.0_dp) then
                   rijn2=rijn*rijn
                   rijn_inv=1.0_dp/rijn
                   rijn_inv2=rijn_inv*rijn_inv
                   rijn_inv4=rijn_inv2*rijn_inv2
                   rijn_inv5=rijn_inv4*rijn_inv
                   damp=1.0_dp/(1.0_dp+exp(-damping_const*((rijn/Svdw0)-1.0_dp)))
                   rsij=rijn/sigij
                   zt1=erf(rsij)-(2.0_dp/sqrt(pi)*rsij*exp(-(rsij**2)))
                   T=0.0_dp
                   do al=1,3
                      do bt=1,3
                         if (al==bt) then
                            zt2=-(3.0_dp*rij(al)*rij(al)-rijn2)*zt1*rijn_inv5
                            zt3=4.0_dp/sqrt(pi)*(rsij**3)*exp(-(rsij**2))*rij(al)*rij(al)*rijn_inv5
                         else
                            zt2=-(3.0_dp*rij(al)*rij(bt))*zt1*rijn_inv5
                            zt3=4.0_dp/sqrt(pi)*(rsij**3)*exp(-(rsij**2))*rij(al)*rij(bt)*rijn_inv5
                         end if
                         T(al,bt)=(zt2+zt3)
                      end do
                   end do
                   TSR=TSR+T*(1.0_dp-damp)

                end if
             end if
          end do
       end do
    end do

  end subroutine

  !> Derivative of the long-range interaction tensor.
  !>
  !> Accumulates the long-range tensor together with derivatives with respect to
  !> pair coordinates and the range-separation radius used by the force routine.
  subroutine tlrijf(ip,jp,kvkl,lattrans,damping_const,Svdw,TLR,dTLRij,dTLRS)
    use constants, only: twopi
    use modinput, only: input

    implicit None

  ! arguments
        !> Fractional coordinates of atom \(i\).
        real(dp), intent(in) :: ip(3)
        !> Fractional coordinates of atom \(j\).
        real(dp), intent(in) :: jp(3)
        !> Fractional k-point vector.
        real(dp), intent(in) :: kvkl(3)
        !> Number of translated unit cells in each lattice direction.
        integer, intent(in) :: lattrans(3)
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> Range-separation radius from screened van der Waals radii.
        real(dp), intent(in) :: Svdw
        !> Long-range interaction tensor.
        complex(dp), intent(out) :: TLR(3,3)
        !> Derivative of the long-range tensor with respect to atom-pair coordinates.
        complex(dp), intent(out) :: dTLRij(3,3,3)
        !> Derivative of the long-range tensor with respect to the range-separation radius.
        complex(dp), intent(out) :: dTLRS(3,3)

  ! local variables
        integer :: i, t1, t2, t3, al, bt, a, b, c
        real(dp) :: kx, ky, kz, ltransn, rijn, kl, zt1
        real(dp) :: rijn2, rijn_inv, rijn_inv2, rijn_inv4, rijn_inv5
        integer :: shft(3)
        real(dp) :: ij(3), rL0(3), ltrans(3), dxyz(3), rij(3)
        real(dp) :: damp, ddamp
        complex(dp) :: T(3,3)
        complex(dp) :: dTlr(3,3,3), dTabc(3,3,3)

    kx=twopi*kvkl(1)
    ky=twopi*kvkl(2)
    kz=twopi*kvkl(3)

    ij=ip-jp
    shft = -nint(ij)
    ij = ij + shft

    rL0=ij(1)*input%structure%crystal%basevect(:,1)+ij(2)*input%structure%crystal%basevect(:,2)&
        &+ij(3)*input%structure%crystal%basevect(:,3)

    TLR=0.0_dp
    T=0.0_dp
    dTLRij=0.0_dp
    dTLRS=0.0_dp
    do t1=-lattrans(1),lattrans(1)
       do t2=-lattrans(2),lattrans(2)
          do t3=-lattrans(3),lattrans(3)
             ltrans=t1*(input%structure%crystal%basevect(:,1))+t2*(input%structure%crystal%basevect(:,2))&
             &+t3*(input%structure%crystal%basevect(:,3))
             ltransn=sqrt(dot_product(ltrans,ltrans))
             dxyz(1)=ij(1)+t1
             dxyz(2)=ij(2)+t2
             dxyz(3)=ij(3)+t3
             if (ltransn<=rsrcutoff) then
                rij=rL0+ltrans
                rijn=sqrt(dot_product(rij,rij))
                if (rijn>0.0_dp) then
                   rijn2=rijn*rijn
                   rijn_inv=1.0_dp/rijn
                   rijn_inv2=rijn_inv*rijn_inv
                   rijn_inv4=rijn_inv2*rijn_inv2
                   rijn_inv5=rijn_inv4*rijn_inv
                   damp=1.0_dp/(1.0_dp+exp(-damping_const*((rijn/Svdw)-1.0_dp)))
                   ddamp=damping_const/Svdw*exp(-damping_const*((rijn/Svdw)-1.0_dp))*damp**2
                   kl=kx*(t1+shft(1))+ky*(t2+shft(2))+kz*(t3+shft(3))
                   T=0.0_dp
                   do al=1,3
                      do bt=1,3
                         if (al==bt) then
                            zt1=(3.0_dp*rij(al)*rij(al)-rijn2)*rijn_inv5
                         else
                            zt1=3.0_dp*rij(al)*rij(bt)*rijn_inv5
                         end if
                         T(al,bt)=zt1
                      end do
                   end do
                   dTlr=0.0_dp
                   do a=1,3
                      do b=1,3
                         do c=1,3
                            if ((a==b).and.(b==c))then
                               dTlr(a,b,c)=derivTlr_a(rij(a),rijn)
                            elseif (a==b) then
                               dTlr(a,b,c)=derivTlr_ab(rij(c),rij(a),rijn)
                            elseif (a==c) then
                               dTlr(a,b,c)=derivTlr_ab(rij(b),rij(a),rijn)
                            elseif (b==c) then
                               dTlr(a,b,c)=derivTlr_ab(rij(a),rij(b),rijn)
                            elseif ((a/=b).and.(b/=c).and.(a/=c)) then
                               dTlr(a,b,c)=derivTlr_abc(rij(a),rij(b),rij(c),rijn)
                            end if
                         end do
                      end do
                   end do
                   dTabc=0.0_dp
                   do i=1,3
                   dTabc(i,:,:)=(dTlr(i,:,:)*(damp)+T(:,:)*ddamp*rij(i)*rijn_inv)*exp(cmplx(0.0_dp, -kl, kind=dp))
                   end do
                   dTLRij=dTLRij+dTabc
                   dTLRS=dTLRS-ddamp*rijn/Svdw*T*exp(cmplx(0.0_dp, -kl, kind=dp))
                   TLR=TLR+T*damp*exp(cmplx(0.0_dp, -kl, kind=dp))
                end if
             end if
          end do
       end do
    end do

  end subroutine

  !> Derivative of the short-range interaction tensor.
  !>
  !> Accumulates the short-range tensor and its pair-coordinate derivatives,
  !> including derivatives of both the tensor and damping function.
  subroutine tsrijf(ip,jp,lattrans,sigij,damping_const,Svdw0,TSR,dTSRij)
    use modmain, only: pi
    use modinput, only: input

    implicit None

  ! arguments
        !> Fractional coordinates of atom \(i\).
        real(dp), intent(in) :: ip(3)
        !> Fractional coordinates of atom \(j\).
        real(dp), intent(in) :: jp(3)
        !> Number of translated unit cells in each lattice direction.
        integer, intent(in) :: lattrans(3)
        !> Gaussian attenuation length.
        real(dp), intent(in) :: sigij
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> Range-separation radius from unscreened atoms-in-molecule van der Waals radii.
        real(dp), intent(in) :: Svdw0
        !> Short-range interaction tensor.
        real(dp), intent(out) :: TSR(3,3)
        !> Derivative of the short-range tensor with respect to atom-pair coordinates.
        real(dp), intent(out) :: dTSRij(3,3,3)
  ! local variables
        integer :: i, t1, t2, t3, al, bt, a, b, c
        real(dp) :: ij(3), rL0(3), dxyz(3), ltrans(3), rij(3), T(3,3)
        real(dp) :: dTsr(3,3,3), dTabc(3,3,3)
        real(dp) :: ltransn, rijn, rsij, zt1, zt2, zt3
        real(dp) :: rijn2, rijn_inv, rijn_inv2, rijn_inv4, rijn_inv5
        real(dp) :: damp,ddamp

    ij=ip-jp
    ij = ij - nint(ij)

    rL0=ij(1)*input%structure%crystal%basevect(:,1)+ij(2)*input%structure%crystal%basevect(:,2)&
        &+ij(3)*input%structure%crystal%basevect(:,3)

    TSR=0.0_dp
    T=0.0_dp
    dTSRij=0.0_dp
    do t1=-lattrans(1),lattrans(1)
       do t2=-lattrans(2),lattrans(2)
          do t3=-lattrans(3),lattrans(3)
             ltrans=t1*(input%structure%crystal%basevect(:,1))+t2*(input%structure%crystal%basevect(:,2))&
             &+t3*(input%structure%crystal%basevect(:,3))
             ltransn=sqrt(dot_product(ltrans,ltrans))
             dxyz(1)=ij(1)+t1
             dxyz(2)=ij(2)+t2
             dxyz(3)=ij(3)+t3
             if (ltransn<=rsrcutoff) then
                rij=rL0+ltrans
                rijn=sqrt(dot_product(rij,rij))
                if (rijn>0.0_dp) then
                   rijn2=rijn*rijn
                   rijn_inv=1.0_dp/rijn
                   rijn_inv2=rijn_inv*rijn_inv
                   rijn_inv4=rijn_inv2*rijn_inv2
                   rijn_inv5=rijn_inv4*rijn_inv
                   damp=1.0_dp/(1.0_dp+exp(-damping_const*((rijn/Svdw0)-1.0_dp)))
                   ddamp=damping_const/Svdw0*exp(-damping_const*((rijn/Svdw0)-1.0_dp))*damp**2
                   rsij=rijn/sigij
                   zt1=erf(rsij)-(2.0_dp/sqrt(pi)*rsij*exp(-(rsij**2)))
                   T=0.0_dp
                   do al=1,3
                      do bt=1,3
                         if (al==bt) then
                            zt2=-(3.0_dp*rij(al)*rij(al)-rijn2)*zt1*rijn_inv5
                            zt3=4.0_dp/sqrt(pi)*(rsij**3)*exp(-(rsij**2))*rij(al)*rij(al)*rijn_inv5
                         else
                            zt2=-(3.0_dp*rij(al)*rij(bt))*zt1*rijn_inv5
                            zt3=4.0_dp/sqrt(pi)*(rsij**3)*exp(-(rsij**2))*rij(al)*rij(bt)*rijn_inv5
                         end if
                         T(al,bt)=(zt2+zt3)
                      end do
                   end do
                   dTsr=0.0_dp
                   do a=1,3
                      do b=1,3
                         do c=1,3
                            if ((a==b).and.(b==c))then
                               dTsr(a,b,c)=derivTsr_a(rij(a),rijn,sigij)
                            elseif (a==b) then
                               dTsr(a,b,c)=derivTsr_ab(rij(c),rij(a),rijn,sigij)
                            elseif (a==c) then
                               dTsr(a,b,c)=derivTsr_ab(rij(b),rij(a),rijn,sigij)
                            elseif (b==c) then
                               dTsr(a,b,c)=derivTsr_ab(rij(a),rij(b),rijn,sigij)
                            elseif ((a/=b).and.(b/=c).and.(a/=c)) then
                               dTsr(a,b,c)=derivTsr_abc(rij(a),rij(b),rij(c),rijn,sigij)
                            end if
                         end do
                      end do
                   end do
                   dTabc=0.0_dp
                   do i=1,3
                      dTabc(i,:,:)=dTsr(i,:,:)*(1.0_dp-damp)-T(:,:)*ddamp*rij(i)*rijn_inv
                   end do
                   dTSRij=dTSRij+dTabc
                   TSR=TSR+T*(1.0_dp-damp)
                end if
             end if
          end do
       end do
    end do

  end subroutine

  !> Compute self-consistently screened polarizabilities using a TS
  !> atoms-in-molecule approach and screened van der Waals radii.
  !>
  !> Hirshfeld volume ratios first scale the free-atom TS parameters. For every
  !> quadrature frequency, the routine constructs and inverts the short-range
  !> response matrix, then takes each atomic block trace to obtain \(\bar\alpha_i\).
  subroutine shortMBD(lattrans, beta, damping_const, alphaiso, Rvdw)
    use modmain, only: sprmax, nspecies, natoms, atposc, spzn, idxas, pi, natmtot
    use modinput, only: input
    use modmpi, only: mpiglobal
    use errors_warnings, only: terminate_if_false
    use TS_vdW_module, only: get_free_atom_vdw_param, atoms_in_ambit, list_of_positions_hirshfeld, list_of_species_hirshfeld, &
                             num_of_atoms_in_sphere_hirshfeld, integrand_numerator, integrand_denominator, current_species,   &
                             current_atom, sph_int

    implicit None

  ! arguments
        !> Number of translated unit cells in each lattice direction.
        integer, intent(in) :: lattrans(3)
        !> Range-separation scaling parameter.
        real(dp), intent(in) :: beta
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> Screened dynamic atoms-in-molecule polarizabilities.
        real(dp), intent(out) :: alphaiso(omegadim,natmtot)
        !> Screened atom-resolved van der Waals radii.
        real(dp), intent(out) :: Rvdw(natmtot)
  ! local variables
        integer :: is, ia, iat, io, i, j, id, jd
        integer :: info
        integer :: nsph, nr
        real(dp) :: C6_free(nspecies), alpha_free_is(nspecies), R0_free(nspecies)
        real(dp) :: R0_eff(natmtot), C6_TS(natmtot), alpha_TS(natmtot)
        real(dp) :: max_sprmax
        real(dp) :: V_eff, V_free, V_ratio
        real(dp) :: alpha_free(natmtot)
        real(dp) :: omegan, sigij, Svdw0
        real(dp) :: omegap(natmtot)
        real(dp) :: alphap(natmtot)
        real(dp) :: sig(natmtot)
        real(dp) :: TSR(3,3)
        real(dp) :: Bij(3*natmtot,3*natmtot)
        real(dp) :: Aij(3*natmtot,3*natmtot)
        real(dp) :: A(3,3)
        real(dp), allocatable :: work(:)
        integer, allocatable :: ipiv(:)

    do is = 1, nspecies
       call get_free_atom_vdw_param(-spzn(is), C6_free(is), alpha_free_is(is), R0_free(is))
    end do

    ! Tkatchenko-Scheffler (TS) atomic polarizability.
    max_sprmax = maxval(sprmax)
    nsph=input%groundstate%TSvdWparameters%nsph
    nr=input%groundstate%TSvdWparameters%nr
    do current_species = 1,nspecies
      do current_atom = 1,natoms(current_species)
        call atoms_in_ambit(max_sprmax + sprmax(current_species), atposc(:, current_atom, current_species), list_of_positions_hirshfeld, list_of_species_hirshfeld, num_of_atoms_in_sphere_hirshfeld)
        V_eff = sph_int(atposc(:, current_atom, current_species), R0_free(current_species), nsph,nr,integrand_numerator)
        V_free = sph_int((/ 0.0_dp, 0.0_dp, 0.0_dp /), R0_free(current_species), 1, 80, integrand_denominator)
        V_ratio = V_eff/V_free
        R0_eff(idxas(current_atom, current_species)) = V_ratio**(1.0_dp/3.0_dp)*R0_free(current_species)
        alpha_TS(idxas(current_atom, current_species)) = V_ratio*alpha_free_is(current_species)
        C6_TS(idxas(current_atom, current_species)) = C6_free(current_species)
        alpha_free(idxas(current_atom, current_species)) = alpha_free_is(current_species)
        end do
    end do

    do iat=1,natmtot
       omegap(iat) = 4.0_dp/3.0_dp*C6_TS(iat)/(alpha_free(iat)**2.0_dp)
    end do

    do io=1,omegadim
       omegan=omg(io)
       Bij=0.0_dp
       do iat=1,natmtot
          alphap(iat)=alpha_TS(iat)/(1.0_dp+(omegan/omegap(iat))**2.0_dp)
          sig(iat)=(sqrt(2.0_dp/pi)*alphap(iat)/3.0_dp)**(1.0_dp/3.0_dp)
       end do

       do i=1,natmtot
          do j=i,natmtot
             sigij=sqrt(sig(i)*sig(i)+sig(j)*sig(j))
             Svdw0=beta*(R0_eff(i)+R0_eff(j))
             call tsrij(xyz(:,i),xyz(:,j),lattrans,sigij,damping_const,Svdw0,TSR)
             do id=1,3
                do jd=1,3
                   if((i==j) .and. (id==jd)) then
                     Bij(3*i-3+id,3*j-3+jd)=1.0_dp/alphap(i)+TSR(id,jd)
                   else
                     Bij(3*i-3+id,3*j-3+jd)=TSR(id,jd)
                     Bij(3*j-3+jd,3*i-3+id)=TSR(id,jd)
                   end if
                end do
             end do
           end do
       end do
       allocate(work(3*natmtot))
       allocate(ipiv(3*natmtot))
       Aij(1:3*natmtot,1:3*natmtot)=Bij(1:3*natmtot,1:3*natmtot)
       call dgetrf(3*natmtot,3*natmtot,Aij,3*natmtot,ipiv,info)
       call terminate_if_false(mpiglobal, info == 0, '(shortMBD): dgetrf failed.')
       call dgetri(3*natmtot,Aij,3*natmtot,ipiv,work,3*natmtot,info)
       call terminate_if_false(mpiglobal, info == 0, '(shortMBD): dgetri failed.')
       deallocate(work,ipiv)

       do i=1,natmtot
          A=0.0_dp
          do j=1,natmtot
             do id=1,3
                do jd=1,3
                   A(id,jd)=A(id,jd)+Aij(3*i-3+id,3*j-3+jd)
                end do
             end do
          end do
          alphaiso(io,i)=(A(1,1)+A(2,2)+A(3,3))/3.0_dp
       end do
    end do

    Rvdw(:)=R0_eff(:)*(alphaiso(omegadim,:)/alpha_TS(:))**(1.0_dp/3.0_dp)

  end subroutine

  !> Compute the MBD@rsSCS dispersion energy.
  !>
  !> Constructs the reciprocal-space long-range interaction matrix, combines it
  !> with the screened polarizabilities, and integrates its eigenvalue response
  !> over imaginary frequency and Brillouin-zone k points.
  subroutine longMBD(beta, damping_const, Edisp)
    use modmain, only: nspecies, natoms, idxas, natmtot, nkpt, vkl, pi, wkpt
    use modinput, only: input
    use modmpi, only: mpi_env_k, distribute_loop
    use errors_warnings, only: terminate_if_false
    use exciting_mpi, only: xmpi_allreduce

    implicit None

  ! arguments
        !> Range-separation scaling parameter.
        real(dp), intent(in) :: beta
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> MBD@rsSCS dispersion energy.
        real(dp), intent(out) :: Edisp
  ! local variables
        integer :: lattrans(3)
        Integer :: is, ia, ias, iat, io, ik, i, j, id, jd, firstk, lastk
        integer :: info
        real(dp) :: alphaiso(omegadim,natmtot)
        real(dp) :: Rvdw(natmtot)
        real(dp) :: Svdw, ed
        real(dp) :: alpha
        real(dp) :: xn(3*natmtot)
        real(dp) :: edisp_buffer(1)
        complex(dp) :: al,bt
        complex(dp) :: TLR(3,3)
        complex(dp) :: ALR(3*natmtot,3*natmtot)
        complex(dp) :: TLRk(3*natmtot,3*natmtot)
        complex(dp) :: AT(3*natmtot,3*natmtot)
        complex(dp) :: evat(3*natmtot)
        complex(dp), allocatable :: vl(:,:)
        complex(dp), allocatable :: vr(:,:)
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)

    rsrcutoff=input%groundstate%MBDparameters%cutoff

    allocate (xyz(3, natmtot), omg(omegadim), omegaweight(omegadim))

    call getlatttrans(lattrans)

    call frequency_grid(omegadim, omg, omegaweight)

    do is=1,nspecies
       do ia=1,natoms(is)
            xyz(:, idxas(ia, is)) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
       end do
    end do

    call shortMBD(lattrans, beta, damping_const, alphaiso, Rvdw)

       Edisp=0.0_dp
#ifdef MPI
    call distribute_loop(mpi_env_k, nkpt, firstk, lastk)
    do ik = firstk, lastk
#else
  do ik=1,nkpt
#endif
       TLRk=0.0_dp
        ed=0.0_dp
       do i=1,natmtot
         do j=i,natmtot
             Svdw=beta*(Rvdw(i)+Rvdw(j))
             call tlrij(xyz(:,i),xyz(:,j),vkl(:, ik),lattrans,damping_const,Svdw,TLR)
             TLRk(3*i-2:3*i, 3*j-2:3*j) = TLR
             if (i.ne.j) then
                TLRk(3*j-2:3*j, 3*i-2:3*i) = conjg(TLR)
             end if
          end do
       end do
       do io=1,omegadim
          ALR=0.0_dp
          do i=1,natmtot
             alpha=alphaiso(io,i)
             ALR(3*i-2,3*i-2)=cmplx(alpha, kind=dp)
             ALR(3*i-1,3*i-1)=cmplx(alpha, kind=dp)
             ALR(3*i,3*i)=cmplx(alpha, kind=dp)
          end do

          AT=-matmul(ALR,TLRk)
          allocate(vl(3*natmtot,3*natmtot),vr(3*natmtot,3*natmtot),work(12*natmtot),rwork(6*natmtot))
          call zgeev('n','n',3*natmtot,AT,3*natmtot,evat,vl,3*natmtot,vr,3*natmtot,work,12*natmtot,rwork,info)
          call terminate_if_false(mpi_env_k, info == 0, '(longMBD): zgeev failed.')
          deallocate(vl,vr,work,rwork)

          if (input%groundstate%MBDparameters%evalremap) then
          ! Eigenvalue remapping (T. Gould et al., J. Chem. Theory Comput. 12, 2016).
             do i=1,3*natmtot
                if (real(evat(i), dp)>=0.0_dp) then
                   xn(i)=evat(i)
                else if (real(evat(i), dp)<0.0_dp) then
                   xn(i)=-(erf(((sqrt(pi)/2.0_dp)*abs(real(evat(i), dp)))**4.0_dp))**0.25_dp
                end if
                   ed=ed-omegaweight(io)*(log(1.0_dp+xn(i))-(xn(i)))
             end do
          else
             do i=1,3*natmtot
                if (real(evat(i), dp)>(-1.0_dp)) then
                    xn(i)=evat(i)
                else
                   call terminate_if_false(mpi_env_k, .false., '(longMBD): polarization catastrophe detected. Set MBDparameters evalremap="true".')
                end if
                   ed=ed-omegaweight(io)*(log(1.0_dp+xn(i)))
             end do
          end if
       end do
       Edisp=Edisp-ed/(2.0_dp*pi)*wkpt(ik)
    end do

#ifdef MPI
    edisp_buffer(1) = Edisp
    call xmpi_allreduce(edisp_buffer, mpi_env_k)
    Edisp = edisp_buffer(1)
#endif
    deallocate(xyz, omg, omegaweight)

  end subroutine

  !> Compute derivatives of self-consistently screened polarizabilities
  !> using a TS atoms-in-molecule approach.
  !>
  !> Differentiates the inverse short-range response matrix to obtain Cartesian
  !> derivatives of each atom-resolved isotropic polarizability.
  subroutine shortfMBD(lattrans, beta, damping_const, alphaiso, Rvdw, dalphaiso)
    use modmain, only: sprmax, nspecies, natoms, atposc, spzn, idxas, pi, natmtot
    use modinput, only: input
    use modmpi, only: mpiglobal
    use errors_warnings, only: terminate_if_false
    use TS_vdW_module, only: get_free_atom_vdw_param, atoms_in_ambit, list_of_positions_hirshfeld, list_of_species_hirshfeld, &
                             num_of_atoms_in_sphere_hirshfeld, integrand_numerator, integrand_denominator, current_species,   &
                             current_atom, sph_int

    implicit None

  ! arguments
        !> Number of translated unit cells in each lattice direction.
        integer, intent(in) :: lattrans(3)
        !> Range-separation scaling parameter.
        real(dp), intent(in) :: beta
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> Screened dynamic atoms-in-molecule polarizabilities.
        real(dp), intent(out) :: alphaiso(omegadim,natmtot)
        !> Screened atom-resolved van der Waals radii.
        real(dp), intent(out) :: Rvdw(natmtot)
        !> Derivative of the screened dynamic polarizabilities.
        real(dp), intent(out) :: dalphaiso(omegadim,natmtot,3*natmtot)
  ! local variables
        integer :: is, ia, iat, io, i, j, k, id, jd, l, m
        integer :: info
        integer :: nsph, nr
        real(dp) :: C6_free(nspecies), alpha_free_is(nspecies), R0_free(nspecies)
        real(dp) :: R0_eff(natmtot), C6_TS(natmtot), alpha_TS(natmtot)
        real(dp) :: max_sprmax
        real(dp) :: V_eff, V_free, V_ratio
        real(dp) :: alpha_free(natmtot)
        real(dp) :: omegan, sigij, Svdw0
        real(dp) :: omegap(natmtot)
        real(dp) :: alphap(natmtot)
        real(dp) :: sig(natmtot)
        real(dp) :: TSR(3,3)
        real(dp) :: dTSRij(3,3,3)
        real(dp) :: Bij(3*natmtot,3*natmtot)
        real(dp) :: Aij(3*natmtot,3*natmtot)
        real(dp) :: dT(3,3*natmtot,3*natmtot)
        real(dp) :: dTSR(3*natmtot,3*natmtot)
        real(dp) :: dAij(3*natmtot,3*natmtot)
        real(dp) :: A(3,3)
        real(dp) :: dA(3,3)
        real(dp), allocatable :: work(:)
        integer, allocatable :: ipiv(:)

    do is = 1, nspecies
       call get_free_atom_vdw_param(-spzn(is), C6_free(is), alpha_free_is(is), R0_free(is))
    end do

    ! Tkatchenko-Scheffler (TS) atomic polarizability.
    max_sprmax = maxval(sprmax)
    nsph=input%groundstate%TSvdWparameters%nsph
    nr=input%groundstate%TSvdWparameters%nr
    do current_species = 1,nspecies
      do current_atom = 1,natoms(current_species)
        call atoms_in_ambit(max_sprmax + sprmax(current_species), atposc(:, current_atom, current_species), list_of_positions_hirshfeld, list_of_species_hirshfeld, num_of_atoms_in_sphere_hirshfeld)
        V_eff = sph_int(atposc(:, current_atom, current_species), R0_free(current_species), nsph,nr,integrand_numerator)
        V_free = sph_int((/ 0.0_dp, 0.0_dp, 0.0_dp /), R0_free(current_species), 1, 80, integrand_denominator)
        V_ratio = V_eff/V_free
        R0_eff(idxas(current_atom, current_species)) = V_ratio**(1.0_dp/3.0_dp) * R0_free(current_species)
        alpha_TS(idxas(current_atom, current_species))= V_ratio*alpha_free_is(current_species)
        C6_TS(idxas(current_atom, current_species)) = C6_free(current_species)
        alpha_free(idxas(current_atom, current_species)) = alpha_free_is(current_species)
        end do
    end do

    do iat=1,natmtot
       omegap(iat) = 4.0_dp/3.0_dp*C6_TS(iat)/(alpha_free(iat)**2.0_dp)
    end do

    dAij=0.0_dp
    do io=1,omegadim
       omegan=omg(io)
       Bij=0.0_dp
       dT=0.0_dp
       do iat=1,natmtot
          alphap(iat)=alpha_TS(iat)/(1.0_dp+(omegan/omegap(iat))**2.0_dp)
          sig(iat)=(sqrt(2.0_dp/pi)*alphap(iat)/3.0_dp)**(1.0_dp/3.0_dp)
       end do

#ifdef USEOMP
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP& SHARED(natmtot,sig,beta,R0_eff,xyz,lattrans,damping_const,Bij,alphap,dT) &
!$OMP& PRIVATE(i,j,id,jd,sigij,Svdw0,TSR,dTSRij)
#endif
       do i=1,natmtot
          do j=i,natmtot
             sigij=sqrt(sig(i)*sig(i)+sig(j)*sig(j))
             Svdw0=beta*(R0_eff(i)+R0_eff(j))
             call tsrijf(xyz(:,i),xyz(:,j),lattrans,sigij,damping_const,Svdw0,TSR,dTSRij)
             do id=1,3
                do jd=1,3
                   if((i==j) .and. (id==jd)) then
                     Bij(3*i-3+id,3*j-3+jd)=1.0_dp/alphap(i)+TSR(id,jd)
                   else
                     Bij(3*i-3+id,3*j-3+jd)=TSR(id,jd)
                     Bij(3*j-3+jd,3*i-3+id)=TSR(id,jd)
                   end if
                end do
             end do
             dT(1:3,3*j-2:3*j,3*i-2)=dTSRij(1,:,:)
             dT(1:3,3*j-2:3*j,3*i-1)=dTSRij(2,:,:)
             dT(1:3,3*j-2:3*j,3*i  )=dTSRij(3,:,:)
             if (i.ne.j) then
             dT(1:3,3*i-2:3*i,3*j-2)=-dTSRij(1,:,:)
             dT(1:3,3*i-2:3*i,3*j-1)=-dTSRij(2,:,:)
             dT(1:3,3*i-2:3*i,3*j  )=-dTSRij(3,:,:)
             end if
           end do
       end do
#ifdef USEOMP
!$OMP END PARALLEL DO
#endif

       allocate(work(3*natmtot))
       allocate(ipiv(3*natmtot))
       Aij(1:3*natmtot,1:3*natmtot)=Bij(1:3*natmtot,1:3*natmtot)
       call dgetrf(3*natmtot,3*natmtot,Aij,3*natmtot,ipiv,info)
       call terminate_if_false(mpiglobal, info == 0, '(shortfMBD): dgetrf failed.')
       call dgetri(3*natmtot,Aij,3*natmtot,ipiv,work,3*natmtot,info)
       call terminate_if_false(mpiglobal, info == 0, '(shortfMBD): dgetri failed.')
       deallocate(work,ipiv)

       do i=1,natmtot
          A=0.0_dp
          do j=1,natmtot
             do id=1,3
                do jd=1,3
                   A(id,jd)=A(id,jd)+Aij(3*i-3+id,3*j-3+jd)
                end do
             end do
          end do
          alphaiso(io,i)=(A(1,1)+A(2,2)+A(3,3))/3.0_dp
       end do

       dalphaiso(io,:,:)=0.0_dp
       do i=1,3*natmtot
          dTSR=0.0_dp
          j=(i-1)/3+1
          do k=1,natmtot
             dTSR(3*j-2:3*j,3*k-2:3*k)=dT(1:3,3*k-2:3*k,i)
             dTSR(3*k-2:3*k,3*j-2:3*j)=dT(1:3,3*k-2:3*k,i)
          end do

          dAij=matmul(Aij,dTSR)
          dAij=matmul(dAij,Aij)

          do l=1,natmtot
             dA=0.0_dp
             do m=1,natmtot
                do id=1,3
                   do jd=1,3
                      dA(id,jd)=dA(id,jd)+dAij(3*l-3+id,3*m-3+jd)
                   end do
                end do
             end do
             dalphaiso(io,l,i)=(dA(1,1)+dA(2,2)+dA(3,3))/3.0_dp
          end do
       end do
    end do

    Rvdw(:)=R0_eff(:)*(alphaiso(omegadim,:)/alpha_TS(:))**(1.0_dp/3.0_dp)

  end subroutine

  !> Compute interatomic MBD@rsSCS forces from the dispersion-energy gradient.
  !>
  !> Combines derivatives of the screened polarizabilities, long-range tensor,
  !> and screening radii, and evaluates the trace of the differentiated response
  !> matrix for every frequency and k point.
  subroutine longfMBD(beta, damping_const, Fdisp)
    use modmain, only: nspecies, natoms, idxas, natmtot, nkpt, vkl, pi, wkpt
    use modinput, only: input
    use modmpi, only: mpi_env_k, distribute_loop
    use errors_warnings, only: terminate_if_false
    use exciting_mpi, only: xmpi_allreduce

    implicit None

  ! arguments
        !> Range-separation scaling parameter.
        real(dp), intent(in) :: beta
        !> Damping parameter used for the Fermi-type damping function.
        real(dp), intent(in) :: damping_const
        !> MBD@rsSCS dispersion forces for all atoms.
        real(dp), intent(out) :: Fdisp(3,natmtot)
  ! local variables
        integer :: lattrans(3)
        Integer :: is, ia, io, ik, i, j, id, jd, l, m, firstk, lastk
        integer :: info
        real(dp) :: alphaiso(omegadim,natmtot)
        real(dp) :: Rvdw(natmtot)
        real(dp) :: dalphaiso(omegadim,natmtot,3*natmtot)
        real(dp) :: Svdw
        real(dp) :: alpha, dSvdwr
        complex(dp) :: TLR(3,3)
        complex(dp) :: dT(3,3*natmtot,3*natmtot)
        complex(dp) :: dTS(3*natmtot,3*natmtot)
        complex(dp) :: dTLRij(3,3,3)
        complex(dp) :: dTLRS(3,3)
        complex(dp) :: ALR(3*natmtot,3*natmtot)
        complex(dp) :: TLRk(3*natmtot,3*natmtot)
        complex(dp) :: ATLR(3*natmtot,3*natmtot)
        complex(dp) :: AT(3*natmtot,3*natmtot)
        complex(dp) :: invAT(3*natmtot,3*natmtot)
        complex(dp) :: dALR(3*natmtot,3*natmtot)
        complex(dp) :: dTLR(3*natmtot,3*natmtot)
        complex(dp) :: Gkw(3*natmtot,3*natmtot)
        complex(dp) :: Fd(3*natmtot,3*natmtot)
        real(dp) :: fforce(3,natmtot)
        complex(dp), allocatable :: work(:)
        integer, allocatable :: ipiv(:)

    rsrcutoff=input%groundstate%MBDparameters%cutoff

    allocate(xyz(3, natmtot), omg(omegadim), omegaweight(omegadim))

    call getlatttrans(lattrans)

    call frequency_grid(omegadim, omg, omegaweight)

    do is=1,nspecies
       do ia=1,natoms(is)
            xyz(:, idxas(ia, is)) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
       end do
    end do

    call shortfMBD(lattrans, beta, damping_const, alphaiso, Rvdw, dalphaiso)

    Fdisp=0.0_dp
#ifdef MPI
    call distribute_loop(mpi_env_k, nkpt, firstk, lastk)
    do ik = firstk, lastk
#else
  do ik=1,nkpt
#endif
       TLRk=0.0_dp
         dT=0.0_dp
        dTS=0.0_dp
      fforce=0.0_dp
         Gkw=0.0_dp
          Fd=0.0_dp
#ifdef USEOMP
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP& SHARED(natmtot,beta,Rvdw,xyz,vkl,ik,lattrans,damping_const,TLRk,dTS,dT) &
!$OMP& PRIVATE(i,j,Svdw,TLR,dTLRij,dTLRS)
#endif
       do i=1,natmtot
          do j=i,natmtot
             Svdw=beta*(Rvdw(i)+Rvdw(j))
             call tlrijf(xyz(:,i),xyz(:,j),vkl(:, ik),lattrans,damping_const,Svdw,TLR,dTLRij,dTLRS)
             TLRk(3*i-2:3*i, 3*j-2:3*j) = TLR
             if (i/=j) then
                TLRk(3*j-2:3*j, 3*i-2:3*i) = conjg(TLR)
             end if
             dTS(3*i-2:3*i,3*j-2)=dTLRS(:,1)
             dTS(3*i-2:3*i,3*j-1)=dTLRS(:,2)
             dTS(3*i-2:3*i,3*j  )=dTLRS(:,3)
             if (i.ne.j) then
                dTLRS=conjg(dTLRS)
                dTS(3*j-2:3*j,3*i-2)=dTLRS(:,1)
                dTS(3*j-2:3*j,3*i-1)=dTLRS(:,2)
                dTS(3*j-2:3*j,3*i  )=dTLRS(:,3)
             end if
             dTLRij(1,:,:)=conjg(dTLRij(1,:,:))
             dTLRij(2,:,:)=conjg(dTLRij(2,:,:))
             dTLRij(3,:,:)=conjg(dTLRij(3,:,:))
  !          dTLRij=conjg(dTLRij)
             dT(1:3,3*j-2:3*j,3*i-2)=dTLRij(1,:,:)
             dT(1:3,3*j-2:3*j,3*i-1)=dTLRij(2,:,:)
             dT(1:3,3*j-2:3*j,3*i  )=dTLRij(3,:,:)
             if (i.ne.j) then
             dT(1:3,3*i-2:3*i,3*j-2)=-dTLRij(1,:,:)
             dT(1:3,3*i-2:3*i,3*j-1)=-dTLRij(2,:,:)
             dT(1:3,3*i-2:3*i,3*j  )=-dTLRij(3,:,:)
             end if
          end do
       end do
#ifdef USEOMP
!$OMP END PARALLEL DO
#endif

       ATLR=0.0_dp
       do i=1,3*natmtot
          ATLR(i,i)=cmplx(1.0_dp, kind=dp)
       end do

       do io=1,(omegadim)
          ALR=0.0_dp
          do i=1,natmtot
             alpha=alphaiso(io,i)
             ALR(3*i-2,3*i-2)=cmplx(alpha, kind=dp)
             ALR(3*i-1,3*i-1)=cmplx(alpha, kind=dp)
             ALR(3*i,3*i)=cmplx(alpha, kind=dp)
          end do

          AT=ATLR-matmul(ALR,TLRk)

          allocate(work(3*natmtot))
          allocate(ipiv(3*natmtot))
          invAT(:,:)=AT(:,:)
          call zgetrf(3*natmtot,3*natmtot,invAT,3*natmtot,ipiv,info)
          call terminate_if_false(mpi_env_k, info == 0, '(longfMBD): zgetrf failed.')
          call zgetri(3*natmtot,invAT,3*natmtot,ipiv,work,3*natmtot,info)
          call terminate_if_false(mpi_env_k, info == 0, '(longfMBD): zgetri failed.')
          deallocate(work,ipiv)

          do i=1,natmtot
             j=3*i-2
             dALR=0.0_dp
             do l=1,natmtot
                dALR(3*l-2,3*l-2)=cmplx(dalphaiso(io,l,j), kind=dp)
                dALR(3*l-1,3*l-1)=cmplx(dalphaiso(io,l,j), kind=dp)
                dALR(3*l,3*l)=cmplx(dalphaiso(io,l,j), kind=dp)
             end do

             dALR=matmul(dALR,TLRk)

             dTLR=0.0_dp
             do l=1,natmtot
                if(i.ge.l) then
                dTLR(3*i-2:3*i,3*l-2:3*l)=dT(1:3,3*l-2:3*l,j)
                dTLR(3*l-2:3*l,3*i-2:3*i)=conjg(dT(1:3,3*l-2:3*l,j))
                else
                dTLR(3*i-2:3*i,3*l-2:3*l)=conjg(dT(1:3,3*l-2:3*l,j))
                dTLR(3*l-2:3*l,3*i-2:3*i)=dT(1:3,3*l-2:3*l,j)
                end if
             end do

             do l=1,natmtot
                do m=l,natmtot
                   dSvdwr=beta/3.0_dp*((Rvdw(l)/alphaiso(omegadim,l)*dalphaiso(omegadim,l,j))&
                   &+(Rvdw(m)/alphaiso(omegadim,m)*dalphaiso(omegadim,m,j)))
                   do id=1,3
                      do jd=1,3
                         dTLR(3*l-3+id,3*m-3+jd)=dTLR(3*l-3+id,3*m-3+jd)+dTS(3*l-3+id,3*m-3+jd)*dSvdwr
                         dTLR(3*m-3+jd,3*l-3+id)=dTLR(3*m-3+jd,3*l-3+id)+dTS(3*m-3+jd,3*l-3+id)*dSvdwr
                      end do
                   end do
                end do
             end do

             dTLR=matmul(ALR,dTLR)
             Gkw=dALR+dTLR
             Fd=matmul(invAT,Gkw)

             do l=1,3*natmtot
                fforce(1,i)=fforce(1,i)+omegaweight(io)*(real(Fd(l,l), dp))
             end do
          end do

          do i=1,natmtot
             j=3*i-1
             dALR=0.0_dp
             do l=1,natmtot
                dALR(3*l-2,3*l-2)=cmplx(dalphaiso(io,l,j), kind=dp)
                dALR(3*l-1,3*l-1)=cmplx(dalphaiso(io,l,j), kind=dp)
                dALR(3*l,3*l)=cmplx(dalphaiso(io,l,j), kind=dp)
             end do

             dALR=matmul(dALR,TLRk)

             dTLR=0.0_dp
             do l=1,natmtot
                if(i>=l) then
                dTLR(3*i-2:3*i,3*l-2:3*l)=dT(1:3,3*l-2:3*l,j)
                dTLR(3*l-2:3*l,3*i-2:3*i)=conjg(dT(1:3,3*l-2:3*l,j))
                else
                dTLR(3*i-2:3*i,3*l-2:3*l)=conjg(dT(1:3,3*l-2:3*l,j))
                dTLR(3*l-2:3*l,3*i-2:3*i)=dT(1:3,3*l-2:3*l,j)
                end if
             end do

             do l=1,natmtot
                do m=l,natmtot
                   dSvdwr=beta/3.0_dp*((Rvdw(l)/alphaiso(omegadim,l)*dalphaiso(omegadim,l,j))&
                   &+(Rvdw(m)/alphaiso(omegadim,m)*dalphaiso(omegadim,m,j)))
                   do id=1,3
                      do jd=1,3
                         dTLR(3*l-3+id,3*m-3+jd)=dTLR(3*l-3+id,3*m-3+jd)+dTS(3*l-3+id,3*m-3+jd)*dSvdwr
                         dTLR(3*m-3+jd,3*l-3+id)=dTLR(3*m-3+jd,3*l-3+id)+dTS(3*m-3+jd,3*l-3+id)*dSvdwr
                      end do
                   end do
                end do
             end do

             dTLR=matmul(ALR,dTLR)
             Gkw=dALR+dTLR
             Fd=matmul(invAT,Gkw)

             do l=1,3*natmtot
                fforce(2,i)=fforce(2,i)+omegaweight(io)*(real(Fd(l,l), dp))
             end do
          end do

          do i=1,natmtot
             j=3*i
             dALR=0.0_dp
             do l=1,natmtot
                dALR(3*l-2,3*l-2)=cmplx(dalphaiso(io,l,j), kind=dp)
                dALR(3*l-1,3*l-1)=cmplx(dalphaiso(io,l,j), kind=dp)
                dALR(3*l,3*l)=cmplx(dalphaiso(io,l,j), kind=dp)
             end do

             dALR=matmul(dALR,TLRk)

             dTLR=0.0_dp
             do l=1,natmtot
                if(i>=l) then
                dTLR(3*i-2:3*i,3*l-2:3*l)=dT(1:3,3*l-2:3*l,j)
                dTLR(3*l-2:3*l,3*i-2:3*i)=conjg(dT(1:3,3*l-2:3*l,j))
                else
                dTLR(3*i-2:3*i,3*l-2:3*l)=conjg(dT(1:3,3*l-2:3*l,j))
                dTLR(3*l-2:3*l,3*i-2:3*i)=dT(1:3,3*l-2:3*l,j)
                end if
             end do


             do l=1,natmtot
                do m=l,natmtot
                   dSvdwr=beta/3.0_dp*((Rvdw(l)/alphaiso(omegadim,l)*dalphaiso(omegadim,l,j))&
                   &+(Rvdw(m)/alphaiso(omegadim,m)*dalphaiso(omegadim,m,j)))
                   do id=1,3
                      do jd=1,3
                         dTLR(3*l-3+id,3*m-3+jd)=dTLR(3*l-3+id,3*m-3+jd)+dTS(3*l-3+id,3*m-3+jd)*dSvdwr
                         dTLR(3*m-3+jd,3*l-3+id)=dTLR(3*m-3+jd,3*l-3+id)+dTS(3*m-3+jd,3*l-3+id)*dSvdwr
                      end do
                   end do
                end do
             end do

             dTLR=matmul(ALR,dTLR)
             Gkw=dALR+dTLR
             Fd=matmul(invAT,Gkw)

             do l=1,3*natmtot
                fforce(3,i)=fforce(3,i)+omegaweight(io)*(real(Fd(l,l), dp))
             end do
          end do
       end do
       Fdisp=Fdisp+fforce/(2.0_dp*pi)*wkpt(ik)
    end do
#ifdef MPI
    call xmpi_allreduce(Fdisp, mpi_env_k)
#endif
    deallocate(xyz, omg, omegaweight)

  end subroutine
end module
