Module TS_vdW_module
  Implicit none
!  Real(8), Parameter :: s6=1d0
!  Real(8), Parameter :: sr6=0.94d0
!  Real(8), Parameter :: damping_const=20d0
!  Real(8), Parameter :: cutoff=95d0
  Integer :: current_atom, current_species
  Integer :: num_of_atoms_in_sphere_hirshfeld
  Real(8), Allocatable :: list_of_positions_hirshfeld(:,:) ! cartesian coordinates of all atoms within a sphere with radius r  
  Integer, Allocatable :: list_of_species_hirshfeld(:)     ! corresponding list of species (index of species as found in modinput)   
  Real(8), Allocatable :: C6ab(:,:), R0_eff_ab(:,:)

Contains

!  Subroutine get_TS_parameters(C6ab, R0_eff_ab)
  Subroutine get_TS_parameters()!C6ab, R0_eff_ab
    Use mod_atoms, Only: sprmax, atposc, nspecies, natoms, spzn, idxas, natmtot
    Use modinput
    Implicit None
!    Real(8), Intent(out) :: C6ab(:,:), R0_eff_ab(:,:)
    Integer :: nsph, nr
    Real(8) :: I_numerator, I_denominator
    Real(8) :: V_ratio
    Real(8) :: max_sprmax
    Real(8) :: C6_free(nspecies), alpha_free_is(nspecies)
    Integer :: is, iat, jat
    Real(8) :: R0_free(nspecies)
    Real(8) :: C6_eff(natmtot), alpha_free_idxas(natmtot), R0_eff(natmtot)
    Real(8) :: e_TS_vdW
    max_sprmax = maxval(sprmax)

    Do is = 1, nspecies
       Call get_free_atom_vdw_param(-spzn(is), C6_free(is), alpha_free_is(is), R0_free(is))
    End Do

    nsph=input%groundstate%TSvdWparameters%nsph
    nr=input%groundstate%TSvdWparameters%nr
!    nsph=590 !possible numbers are: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3740, 3890, 4334, 4802, 5294, 5810
!    nr=120
!!!! quick calc:
!    nsph=146
!    nr=80
!!!!

    Do current_species = 1, nspecies
       Do current_atom = 1,natoms(current_species)
          Call atoms_in_ambit(max_sprmax + sprmax(current_species), atposc(:, current_atom, current_species), list_of_positions_hirshfeld, list_of_species_hirshfeld, num_of_atoms_in_sphere_hirshfeld)
          I_numerator = sph_int(atposc(:, current_atom, current_species), R0_free(current_species), nsph,nr,integrand_numerator)
          I_denominator = sph_int( (/ 0d0, 0d0, 0d0 /), R0_free(current_species), 1, 80, integrand_denominator)
          V_ratio = I_numerator/I_denominator
          C6_eff(idxas(current_atom, current_species)) = V_ratio**2 * C6_free(current_species)
          alpha_free_idxas(idxas(current_atom, current_species)) =  alpha_free_is(current_species)
          R0_eff(idxas(current_atom, current_species)) = V_ratio**(1d0/3) * R0_free(current_species)
       End Do ! current_atom
    End Do ! current_species

    Do iat = 1, natmtot
       Do jat = 1, natmtot
          C6ab(iat, jat) = 2d0 / ( alpha_free_idxas(jat)/(alpha_free_idxas(iat)*C6_eff(jat)) + alpha_free_idxas(iat)/(alpha_free_idxas(jat)*C6_eff(iat)) )
          R0_eff_ab(iat, jat) = R0_eff(iat) + R0_eff(jat)
       End Do
    End Do
  End Subroutine get_TS_parameters

  Function integrand_numerator(vpc,np)
    Use modmain, Only: lmmaxvr, rhoir, rhomt
    Use mod_atoms, Only: atposc, sprmax
    Use modinput
    Implicit None
    Integer, Intent(in) :: np
    Real(8), Intent(in) :: vpc(:,:)
    Real(8) :: integrand_numerator(np)
    Integer :: lmax
    Real(8) :: inv_basevect(3,3)
    Integer :: ip
    Real(8) :: rvec(3), r, r_array(np)
    Real(8) :: vpc_buffer(3,np)
    Integer :: bookkeeping(np)
    Real(8), Allocatable :: vpl(:,:), vpc_reduced(:,:)
    Real(8), Allocatable :: rho(:)
    Integer :: count_in
    lmax=input%groundstate%lmaxvr
    r_array(:) = 0.
    integrand_numerator(:)=0.

    ! consider only sites with nonzero hirshfeld-weight 
    count_in = 0
    Do ip=1,np
       rvec(:)=atposc(:, current_atom, current_species)-vpc(:, ip)
       r=Sqrt(Dot_product(rvec(:),rvec(:)))
       If ( r .Lt. sprmax(current_species) ) Then
          count_in = count_in + 1
          vpc_buffer(:,count_in) = vpc(:, ip)
          r_array(count_in) = r
          bookkeeping(count_in) = ip
       End If
    End Do
    If ( count_in .Eq. 0 ) return
    Allocate(vpl(3,count_in), vpc_reduced(3,count_in))
    Allocate(rho(count_in))
    vpc_reduced(:,:)=vpc_buffer(:,1:count_in)

    ! convert cartesian coordinates to lattice coordinates
    Call r3minv (input%structure%crystal%basevect, inv_basevect)
    vpl(:,:)=matmul(inv_basevect(:,:),vpc_reduced(:,:))

    ! get density
    call rfarray(lmax,lmmaxvr,rhomt,rhoir,count_in,vpl, rho)

    ! modify density according to the form of the desired integrand  
    rho(:)=rho(:)/hirshfeldweightdenominator(vpc_reduced, count_in)

    ! build correct integrand
    Do ip=1,count_in
       integrand_numerator(bookkeeping(ip))=r_array(ip)**3 * nfree(r_array(ip),current_species) * rho(ip)
    End Do
  End Function integrand_numerator

  Function integrand_denominator(vpc,np)
    !    Use mod_atoms, Only: atposc
    Implicit None
    Integer, Intent(in) :: np
    Real(8), Intent(in) :: vpc(:,:)
    Real(8) :: integrand_denominator(np)
    Integer :: ip
    !    Real(8) :: rvec(3), r, r3
    Real(8) :: r, r3

    Do ip=1,np
       !       rvec(:)=atposc(:, current_atom, current_species)-vpc(:, ip)
       !       r=Sqrt(Dot_product(rvec(:),rvec(:)))
       r=Sqrt(Dot_product(vpc(:, ip), vpc(:, ip)))
       r3=r**3
       integrand_denominator(ip)=r3 * nfree(r,current_species)
    End Do
  End Function Integrand_denominator

  Function sph_int(vorigin,rm,nsph,nr,integrand)
    ! single center integration
    ! using Lebedev's quadrature for the spherical part
    ! and Gauss-Chebyshev quadrature for the radial part
    Implicit None
    Interface
       Function integrand(vpc,n)
         Integer, Intent(in) :: n
         Real(8), Intent(in) :: vpc(:,:)
         Real(8) :: integrand(n)
       End Function integrand
    End Interface
    Real(8) :: sph_int
    Real(8), Intent(in) :: vorigin(3)
    Real(8), Intent(in) :: rm   ! "midpoint" of the radial integration interval
    Integer, Intent(in) :: nsph ! number of Lebedev grid points, possible numbers are: 1, 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3740, 3890, 4334, 4802, 5294, 5810
    Integer, Intent(in) :: nr   ! number of radial grid points 
    Integer :: ntot
    Integer :: ir, icount, jcount
    Integer :: idxmin, idxmax, idxmodulo
    Real(8) :: pi, r, rmin, G, wr, f
    Real(8), allocatable :: v_unitsph(:,:), vpc(:,:)
    Real(8), allocatable :: wsph(:), x(:), fp(:)

    pi=atan(1.)*4
    ntot=nsph*nr
    allocate(v_unitsph(3,nsph), vpc(3,ntot))
    allocate(wsph(nsph), x(nr), fp(ntot))

    ! nsph = 1 for spherical averaged integrands
    If (nsph .Eq. 1) Then
       v_unitsph(:,1)=(/ 1, 0, 0 /)
       wsph(:)=4*pi
    Else
       ! generate Lebedev grid and weights
       call leblaik(nsph, v_unitsph(:,:),wsph(:))
       wsph(:)=wsph(:)*4*pi
    End If


    ! generate radial grid for Gauss-Chebyshev quadrature of the second kind
    do ir=1,nr
       x(ir)=cos(dble(ir)/(nr+1)*pi)
       r=rm*(1+x(ir))/(1-x(ir))
       idxmin=(ir-1)*nsph + 1
       idxmax=idxmin + nsph -1
       vpc(:,idxmin:idxmax) = r*v_unitsph(:,:)
    end do

    ! shift grid to preset center of integration 
    rmin=rm*(1+x(nr))/(1-x(nr))
    do jcount=1,3
       if (abs(vorigin(jcount)) .gt. epsilon(rmin)*rmin/100) then !!! 
          do icount=1,ntot
             vpc(jcount,icount) = vpc(jcount,icount) + vorigin(jcount)
          enddo
       endif
    enddo

    ! get function evaluations at grid points
    fp=integrand(vpc,ntot)

    ! Integration by summing weighted function evaluations
    sph_int=0.
    do ir=1,nr
       idxmin=(ir-1)*nsph + 1
       idxmax=idxmin + nsph -1
       idxmodulo=0
       G=0.
       do icount = idxmin, idxmax
          idxmodulo = idxmodulo + 1
          G = G + wsph(idxmodulo)*fp(icount)
       enddo
       wr=pi/(nr+1)*(sin(dble(ir)/(nr+1)*pi))**2
       f=((1+x(ir))/(1-x(ir))**3)**(1.5)*G
       r=rm*(1+x(ir))/(1-x(ir))
       sph_int = sph_int + wr*f
    enddo
    sph_int = sph_int*2*rm**3
  End Function sph_int

  Function hirshfeldweightdenominator(vpc, np)
    Use mod_atoms, Only: sprmax
    Implicit None
    Real(8), Intent(in) :: vpc(:,:)
    Integer, Intent(in) :: np
    Real(8) :: hirshfeldweightdenominator(np)
    Integer :: ip, icount
    Real(8) :: rvec(3), r, sprmax_
    hirshfeldweightdenominator(:)=0d0
    sprmax_ = sprmax(current_species)
    Do ip=1,np !MPI ???
       Do  icount = 1, num_of_atoms_in_sphere_hirshfeld 
          rvec(:)=list_of_positions_hirshfeld(:,icount) - vpc(:, ip)
          r=Sqrt(Dot_product(rvec(:),rvec(:)))
          If (r .Gt. sprmax_) Cycle
          hirshfeldweightdenominator(ip)=hirshfeldweightdenominator(ip) + nfree(r,list_of_species_hirshfeld(icount))
       End Do ! icount
    End Do ! ip
  End Function Hirshfeldweightdenominator


  Function nfree(r,is)
    Use modinput
    Use mod_atoms, Only: spnr, spr, sprho, sprmax
    Implicit None
    Real(8), Intent(in) :: r  ! distance
    Integer, Intent(in) :: is ! species
    Real(8) :: nfree
    Integer :: np2, ir, ir0, i, j
    Real(8) :: r_
    Real (8), Allocatable :: ya (:), c (:)
    Real (8) :: polynom
    External polynom
    Allocate (ya(input%groundstate%nprad),&
         c(input%groundstate%nprad))
    If (r .Gt. sprmax(is)) Then
       nfree=0.
    Else
       np2 = input%groundstate%nprad / 2
       Do ir = 1, spnr(is)
          If (spr(ir, is) .Ge. r) Then
             If (ir .Le. np2) Then
                ir0 = 1
             Else If (ir .Gt. spnr(is)-np2) Then !or just nfree=0. ?
                ir0 = spnr(is) - &
                     &     input%groundstate%nprad + 1
             Else
                ir0 = ir - np2
             End If
             r_ = Max (r, spr(1, is))
             Do j = 1, input%groundstate%nprad
                i = ir0 + j - 1
                ya(j) = sprho(i,is)
             End Do
             nfree = polynom(0, input%groundstate%nprad, spr(ir0, is), ya, c, r_)
             Exit
          End If
       End Do
    End If
  End Function nfree

  Subroutine atoms_in_ambit(r_max, center_c, list_of_positions, list_of_species, num_of_atoms_in_sphere)
    Use vdw_general_routines, Only: getlatticerepetition
    Use mod_atoms, Only: natmtot, nspecies, natoms, atposc
    Use modinput
    Implicit none
    Real(8), Intent(in) :: r_max                                ! radius (Bohr)
    Real(8), Intent(in) :: center_c(3)                          ! cartesian coordinates of center
    Real(8), Allocatable, Intent(out) :: list_of_positions(:,:) ! cartesian coordinates of all atoms within a sphere with radius r
    Integer, Allocatable, Intent(out) :: list_of_species(:)     ! corresponding list of species (index of species as found in modinput)
    Integer, Intent(out) :: num_of_atoms_in_sphere
    Real(8), Allocatable :: list_of_positions_extended(:,:)
    Integer, Allocatable :: list_of_species_extended(:)
    Integer :: latticerepetition(3)
    Integer :: max_num_of_atoms
    Integer :: icount
    Integer :: is, ia, tau_a, tau_b, tau_c
    Integer :: center_l_integer(3)
    Real(8) :: tau_coeff(3)
    Real(8) :: inv_basevect(3,3)
    Real(8) :: center_l(3), center_c_equivalent(3), tau(3), center_c_integer(3)
    Real(8) :: rvec(3), r
    Call getlatticerepetition(latticerepetition,r_max)
    max_num_of_atoms = natmtot
    Do icount=1,3
       max_num_of_atoms = max_num_of_atoms * (2*latticerepetition(icount)+1)
    End Do

    Allocate(list_of_positions_extended(3,max_num_of_atoms), &
         list_of_species_extended(max_num_of_atoms))

    ! convert cartesian coordinates to lattice coordinates
    Call r3minv (input%structure%crystal%basevect, inv_basevect)
    Call r3mv (inv_basevect(:,:), center_c(:), center_l)

    ! take equivalent site in the "first" unit cell
    Call r3frac (input%structure%epslat, center_l, center_l_integer)

    ! convert back to cartesian coordinates
    Call r3mv (input%structure%crystal%basevect, center_l, center_c_equivalent)
    Call r3mv (input%structure%crystal%basevect, dble(center_l_integer), center_c_integer)

    icount = 0      
    Do is=1,nspecies
       Do ia=1,natoms(is)
          Do tau_a = -latticerepetition(1),latticerepetition(1)
             Do tau_b = -latticerepetition(2),latticerepetition(2)
                Do tau_c = -latticerepetition(3),latticerepetition(3)
                   tau_coeff = dble((/tau_a, tau_b, tau_c/))
                   Call r3mv(input%structure%crystal%basevect,tau_coeff,tau)
                   rvec(:)=atposc(:, ia, is) + tau(:) - center_c_equivalent(:)
                   r=Sqrt(Dot_product(rvec(:),rvec(:)))
                   If (r .Le. r_max) Then
                      icount = icount + 1
                      list_of_positions_extended(:,icount) = atposc(:, ia, is) + tau(:) + center_c_integer(:)
                      list_of_species_extended(icount) = is
                   End If
                End Do ! tau_c
             End Do ! tau_b
          End Do ! tau_a
       End Do ! ia
    End Do ! is
    num_of_atoms_in_sphere = icount
    Allocate(list_of_positions(3,num_of_atoms_in_sphere), &
         list_of_species(num_of_atoms_in_sphere))
    list_of_positions(:,:) = list_of_positions_extended(:,1:num_of_atoms_in_sphere)
    list_of_species(:) = list_of_species_extended(1:num_of_atoms_in_sphere)
  End Subroutine atoms_in_ambit

  subroutine get_free_atom_vdw_param(nucleus, C6, alpha, R0)
    implicit none
    real*8 :: nucleus, C6,alpha, R0

    if (nucleus.eq.1.) then ! H
       alpha=4.500000
       C6=6.500000
       R0=3.100000

    else if (nucleus.eq.2.) then ! He
       alpha=1.380000
       C6=1.460000
       R0=2.650000

    else if (nucleus.eq.3.) then ! Li
       alpha=164.200000
       C6=1387.000000
       R0=4.160000

    else if (nucleus.eq.4.) then ! Be
       alpha=38.000000
       C6=214.000000
       R0=4.170000

    else if (nucleus.eq.5.) then ! B
       alpha=21.000000
       C6=99.500000
       R0=3.890000

    else if (nucleus.eq.6.) then ! C
       alpha=12.000000
       C6=46.600000
       R0=3.590000

    else if (nucleus.eq.7.) then ! N
       alpha=7.400000
       C6=24.200000
       R0=3.340000

    else if (nucleus.eq.8.) then ! O
       alpha=5.400000
       C6=15.600000
       R0=3.190000

    else if (nucleus.eq.9.) then ! F
       alpha=3.800000
       C6=9.520000
       R0=3.040000

    else if (nucleus.eq.10.) then ! Ne
       alpha=2.670000
       C6=6.380000
       R0=2.910000

    else if (nucleus.eq.11.) then ! Na
       alpha=162.700000
       C6=1556.000000
       R0=3.730000

    else if (nucleus.eq.12.) then ! Mg
       alpha=71.000000
       C6=627.000000
       R0=4.270000

    else if (nucleus.eq.13.) then ! Al
       alpha=60.000000
       C6=528.000000
       R0=4.330000

    else if (nucleus.eq.14.) then ! Si
       alpha=37.000000
       C6=305.000000
       R0=4.200000

    else if (nucleus.eq.15.) then ! P
       alpha=25.000000
       C6=185.000000
       R0=4.010000

    else if (nucleus.eq.16.) then ! S
       alpha=19.600000
       C6=134.000000
       R0=3.860000

    else if (nucleus.eq.17.) then ! Cl
       alpha=15.000000
       C6=94.600000
       R0=3.710000

    else if (nucleus.eq.18.) then ! Ar
       alpha=11.100000
       C6=64.300000
       R0=3.550000

    else if (nucleus.eq.19.) then ! K
       alpha=292.900000
       C6=3897.000000
       R0=3.710000

    else if (nucleus.eq.20.) then ! Ca
       alpha=160.000000
       C6=2221.000000
       R0=4.650000

    else if (nucleus.eq.21.) then ! Sc
       alpha=120.000000
       C6=1383.000000
       R0=4.590000

    else if (nucleus.eq.22.) then ! Ti
       alpha=98.000000
       C6=1044.000000
       R0=4.510000

    else if (nucleus.eq.23.) then ! V
       alpha=84.000000
       C6=832.000000
       R0=4.440000

    else if (nucleus.eq.24.) then ! Cr
       alpha=78.000000
       C6=602.000000
       R0=3.990000

    else if (nucleus.eq.25.) then ! Mn
       alpha=63.000000
       C6=552.000000
       R0=3.970000

    else if (nucleus.eq.26.) then ! Fe
       alpha=56.000000
       C6=482.000000
       R0=4.230000

    else if (nucleus.eq.27.) then ! Co
       alpha=50.000000
       C6=408.000000
       R0=4.180000

    else if (nucleus.eq.28.) then ! Ni
       alpha=48.000000
       C6=373.000000
       R0=3.820000

    else if (nucleus.eq.29.) then ! Cu
       alpha=42.000000
       C6=253.000000
       R0=3.760000

    else if (nucleus.eq.30.) then ! Zn
       alpha=40.000000
       C6=284.000000
       R0=4.020000

    else if (nucleus.eq.31.) then ! Ga
       alpha=60.000000
       C6=498.000000
       R0=4.190000

    else if (nucleus.eq.32.) then ! Ge
       alpha=41.000000
       C6=354.000000
       R0=4.200000

    else if (nucleus.eq.33.) then ! As
       alpha=29.000000
       C6=246.000000
       R0=4.110000

    else if (nucleus.eq.34.) then ! Se
       alpha=25.000000
       C6=210.000000
       R0=4.040000

    else if (nucleus.eq.35.) then ! Br
       alpha=20.000000
       C6=162.000000
       R0=3.930000

    else if (nucleus.eq.36.) then ! Kr
       alpha=16.800000
       C6=129.600000
       R0=3.820000

    else if (nucleus.eq.37.) then ! Rb
       alpha=319.200000
       C6=4691.000000
       R0=3.720000

    else if (nucleus.eq.38.) then ! Sr
       alpha=199.000000
       C6=3170.000000
       R0=4.540000

    else if (nucleus.eq.39.) then ! Y
       alpha=126.7370
       C6=1968.580
       R0=4.81510

    else if (nucleus.eq.40.) then ! Zr
       alpha=119.97
       C6=1677.91
       R0=4.53

    else if (nucleus.eq.41.) then ! Nb
       alpha=101.603
       C6=1263.61
       R0=4.2365

    else if (nucleus.eq.42.) then ! Mo
       alpha=88.4225785
       C6=1028.73
       R0=4.099

    else if (nucleus.eq.43.) then ! Tc
       alpha=80.083
       C6=1390.87
       R0=4.076

    else if (nucleus.eq.44.) then ! Ru
       alpha=65.8950
       C6=609.754
       R0=3.99530

    else if (nucleus.eq.45.) then ! Rh
       alpha=56.1
       C6=469.0
       R0=3.95

    else if (nucleus.eq.46.) then ! Pd
       alpha=23.680000
       C6=157.500000
       R0=3.66000

    else if (nucleus.eq.47.) then ! Ag
       alpha=50.600000
       C6=339.000000
       R0=3.820000

    else if (nucleus.eq.48.) then ! Cd
       alpha=39.7
       C6=452.0
       R0=3.99

    else if (nucleus.eq.49.) then ! In
       alpha=70.22000
       C6=707.046000
       R0=4.23198000

    else if (nucleus.eq.50.) then ! Sn
       alpha=55.9500
       C6=587.41700
       R0=4.303000

    else if (nucleus.eq.51.) then ! Sb
       alpha=43.671970 
       C6=459.322
       R0=4.2760

    else if (nucleus.eq.52.) then ! Te
       alpha=37.65
       C6=396.0
       R0=4.22

    else if (nucleus.eq.53.) then ! I
       alpha=35.000000
       C6=385.000000
       R0=4.170000

    else if (nucleus.eq.54.) then ! Xe
       alpha=27.300000
       C6=285.900000
       R0=4.080000

    else if (nucleus.eq.55.) then ! Cs
       alpha=427.12
       C6=6582.08
       R0=3.78

    else if (nucleus.eq.56.) then ! Ba
       alpha=275.0
       C6=5727.0
       R0=4.77

    else if (nucleus.eq.72.) then ! Hf
       alpha=99.52
       C6=1274.8
       R0=4.21

    else if (nucleus.eq.73.) then ! Ta
       alpha=82.53
       C6=1019.92
       R0=4.15

    else if (nucleus.eq.74.) then ! W
       alpha=71.041
       C6=847.93
       R0=4.08

    else if (nucleus.eq.75.) then ! Re
       alpha=63.04
       C6=710.2
       R0=4.02

    else if (nucleus.eq.76.) then ! Os
       alpha=55.055
       C6=596.67
       R0=3.84


    else if (nucleus.eq.77.) then ! Ir
       alpha=42.51
       C6=359.1
       R0=4.00

    else if (nucleus.eq.78.) then ! Pt
       alpha=39.68
       C6=347.1
       R0=3.92

    else if (nucleus.eq.79.) then ! Au
       alpha=36.5
       C6=298.0
       R0=3.86

    else if (nucleus.eq.80.) then ! Hg
       alpha=33.9
       C6=392.0
       R0=3.98

    else if (nucleus.eq.81.) then ! Tl
       alpha=69.92
       C6=717.44
       R0=3.91

    else if (nucleus.eq.82.) then ! Pb
       alpha=61.8
       C6=697.0
       R0=4.31

    else if (nucleus.eq.83.) then ! Bi
       alpha=49.02
       C6=571.0
       R0=4.32

    else if (nucleus.eq.84.) then ! Po
       alpha=45.013
       C6=530.92
       R0=4.097

    else if (nucleus.eq.85.) then ! At
       alpha=38.93
       C6=457.53
       R0=4.07

    else if (nucleus.eq.86.) then ! Rn
       alpha=33.54
       C6=390.63
       R0=4.23

    else
       ! case default
       C6=1.0 
       alpha=1.0
       R0=1.0
       write(*,*) '*************** WARNING***************: VdW parameters not defined for atom: ', & 
            nucleus

       stop

    end if

  end subroutine get_free_atom_vdw_param

End Module TS_vdW_module
