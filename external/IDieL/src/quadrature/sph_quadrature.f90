! Copyright (C) 2020-2024 GreenX library
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

!> @file
!> Contains elements to compute spherical quadratures using a Lebedev grid of order 131
!> For further information on how the points are generated see the following references:
!>   
!>   [1] V.I. Lebedev, \(\textit{"Quadratures on a sphere"}\), 
!>       USSR Computational Mathematics and Mathematical Physics \(\textbf{16}\), 10-24 (1976)
!>
!>   And for the weights:
!> 
!>   [2] V.I. Lebedev, and D.N. Laikov \(\textit{"A quadrature formula for the sphere of the 131st algebraic order of accuracy"}\)
!>       Doklady Mathematics \(\textbf{59}\), 477-481 (1999)
!>


!> This module contains subroutines to compute meshes and weights for integrals over a sphere using Lebedev quadrature
module idiel_sph_quadrature
    
    use idiel_constants,   only: i32, aip, pi, twopi

    implicit none
    
    private 
    public :: compute_angular_mesh_lebedev_131, & 
              compute_angular_mesh_lebedev_41, compute_angular_mesh_lebedev_21

contains

    !> Spherical to cartesian coordinates
    !> @param[in] ang - angular mesh (theta,phi)
    !> @param[in] r - the radius
    !> @result    xyz - the Cartesian coordinates 
    pure function spherical_2_cartesian(ang,r) result(xyz)
        !> The angular mesh (theta,phi)
        real(aip), intent(in) :: ang(:,:)
        !> The spherical radius
        real(aip), intent(in) :: r
        !> The Cartesian mesh
        real(aip), allocatable :: xyz(:,:)

        allocate(xyz(size(ang,1),3))

        xyz(:,1) = r * sin(ang(:,1)) * cos(ang(:,2))
        xyz(:,2) = r * sin(ang(:,1)) * sin(ang(:,2))
        xyz(:,3) = r * cos(ang(:,1))

    end function spherical_2_cartesian

    !> Spherical Cartesian to angular mesh
    !> @param[in] xyz - angular cartesian mesh (x,y,z)
    !> @param[in] r - the radius
    !> @result    ang - the angular coordinates (theta, phi)
    pure function cartesian_2_spherical(xyz,r) result(ang)
        !> The Cartesian mesh (x,y,z)
        real(aip), intent(in) :: xyz(:,:)
        !> The spherical radius
        real(aip), intent(in) :: r
        !> The angular mesh
        real(aip), allocatable :: ang(:,:)

        allocate(ang(size(xyz,1),2), source=0.0_aip)

        ang(:,1) = acos( xyz(:,3) / r )
        where(hypot(xyz(:,1),xyz(:,2)) .gt. 1.0e-12) ang(:,2) = sign(acos( xyz(:,1) / hypot(xyz(:,1),xyz(:,2)) ), xyz(:,2))

    end function cartesian_2_spherical


    !> Compute an angular mesh using and weights using a Lebedev quadrature of order 21
    !> @param[out]  ang  - the angles (theta,phi)
    !> @param[out]  w    - weights
    !> @param[out]  xyz  - the mesh in cartesian coordinates for unitary radius
    subroutine compute_angular_mesh_lebedev_21(ang, w, xyz)

        real(aip), allocatable, intent(out) :: ang(:,:)
        real(aip), allocatable, intent(out) :: w(:)
        real(aip), allocatable, intent(out) :: xyz(:,:)

        ! Number of points
        integer(i32) :: npoints = 170_i32
        integer(i32) :: n = 1_i32

        ! Allocate the space
        if (allocated(xyz)) deallocate(xyz)
        allocate(xyz(npoints,3))
        if (allocated(w)) deallocate(w)
        allocate(w(npoints))

        ! Build the Lebedev mesh using the Lebedev-Laikov values
        ! and octahedral points
        call generate_a1(xyz, w, n, 0.5544842902037365e-2_aip)
        call generate_a2(xyz, w, n, 0.6071332770670752e-2_aip)
        call generate_a3(xyz, w, n, 0.6383674773515093e-2_aip)
        call generate_bk(xyz, w, 0.2551252621114134_aip, n, 0.5183387587747790e-2_aip)
        call generate_bk(xyz, w, 0.6743601460362766_aip, n, 0.6317929009813725e-2_aip)
        call generate_bk(xyz, w, 0.4318910696719410_aip, n, 0.6201670006589077e-2_aip)
        call generate_ck(xyz, w, 0.2613931360335988_aip, n, 0.5477143385137348e-2_aip)
        call generate_dk(xyz, w, 0.4990453161796037_aip, 0.1446630744325115_aip, n, 0.5968383987681156e-2_aip)

        ! Obtain the spherical angular coordinates
        ang = cartesian_2_spherical(xyz, 1.0_aip)

        ! We finally include de 4*pi term in the weight
        w(:) = 4.0_aip * pi * w(:)

    end subroutine compute_angular_mesh_lebedev_21

    !> Compute an angular mesh using and weights using a Lebedev quadrature of order 41
    !> @param[out]  ang  - the angles (theta,phi)
    !> @param[out]  w    - weights
    !> @param[out]  xyz  - the mesh in cartesian coordinates for unitary radius
    subroutine compute_angular_mesh_lebedev_41(ang, w, xyz)

        real(aip), allocatable, intent(out) :: ang(:,:)
        real(aip), allocatable, intent(out) :: w(:)
        real(aip), allocatable, intent(out) :: xyz(:,:)

        ! Number of points
        integer(i32) :: npoints = 590_i32
        integer(i32) :: n = 1_i32

        ! Allocate the space
        if (allocated(xyz)) deallocate(xyz)
        allocate(xyz(npoints,3))
        if (allocated(w)) deallocate(w)
        allocate(w(npoints))

        ! Build the Lebedev mesh using the Lebedev-Laikov values
        ! and octahedral points
        call generate_a1(xyz, w, n, 0.3095121295306187e-3_aip)
        call generate_a3(xyz, w, n, 0.1852379698597489e-2_aip)
        call generate_bk(xyz, w, 0.7040954938227469_aip, n, 0.1871790639277744e-2_aip)
        call generate_bk(xyz, w, 0.6807744066455243_aip, n, 0.1858812585438317e-2_aip)
        call generate_bk(xyz, w, 0.6372546939258752_aip, n, 0.1852028828296213e-2_aip)
        call generate_bk(xyz, w, 0.5044419707800358_aip, n, 0.1846715956151242e-2_aip)
        call generate_bk(xyz, w, 0.4215761784010967_aip, n, 0.1818471778162769e-2_aip)
        call generate_bk(xyz, w, 0.3317920736472123_aip, n, 0.1749564657281154e-2_aip)
        call generate_bk(xyz, w, 0.2384736701421887_aip, n, 0.1617210647254411e-2_aip)
        call generate_bk(xyz, w, 0.1459036449157763_aip, n, 0.1384737234851692e-2_aip)
        call generate_bk(xyz, w, 0.6095034115507196e-1_aip, n, 0.9764331165051050e-3_aip)
        call generate_ck(xyz, w, 0.6116843442009876_aip, n, 0.1857161196774078e-2_aip)
        call generate_ck(xyz, w, 0.3964755348199858_aip, n, 0.1705153996395864e-2_aip)
        call generate_ck(xyz, w, 0.1724782009907724_aip, n, 0.1300321685886048e-2_aip)
        call generate_dk(xyz, w, 0.5610263808622060_aip,0.3518280927733519_aip, n, 0.1842866472905286e-2_aip)
        call generate_dk(xyz, w, 0.4742392842551980_aip,0.2634716655937950_aip, n, 0.1802658934377451e-2_aip)
        call generate_dk(xyz, w, 0.5984126497885380_aip,0.1816640840360209_aip, n, 0.1849830560443660e-2_aip)
        call generate_dk(xyz, w, 0.3791035407695563_aip,0.1720795225656878_aip, n, 0.1713904507106709e-2_aip)
        call generate_dk(xyz, w, 0.2778673190586244_aip,0.8213021581932511e-1_aip, n, 0.1555213603396808e-2_aip)
        call generate_dk(xyz, w, 0.5033564271075117_aip,0.8999205842074875e-1_aip, n, 0.1802239128008525e-2_aip)

        ! Obtain the spherical angular coordinates
        ang = cartesian_2_spherical(xyz, 1.0_aip)

        ! We finally include de 4*pi term in the weight
        w(:) = 4.0_aip * pi * w(:)

    end subroutine compute_angular_mesh_lebedev_41

    !> Compute an angular mesh using and weights using a Lebedev quadrature of order 131
    !> @param[out]  ang  - the angles (theta,phi)
    !> @param[out]  w    - weights
    !> @param[out]  xyz  - the mesh in cartesian coordinates for unitary radius
    subroutine compute_angular_mesh_lebedev_131(ang, w, xyz)

        real(aip), allocatable, intent(out) :: ang(:,:)
        real(aip), allocatable, intent(out) :: w(:)
        real(aip), allocatable, intent(out) :: xyz(:,:)

        ! Number of points 
        integer(i32) :: npoints = 5810_i32
        integer(i32) :: n = 1_i32

        ! Allocate the space
        if (allocated(xyz)) deallocate(xyz)
        allocate(xyz(npoints,3))
        if (allocated(w)) deallocate(w)
        allocate(w(npoints))

        ! Build the Lebedev mesh using the Lebedev-Laikov values
        ! and octahedral points
        call generate_a1(xyz, w, n, 0.9735347946175486e-5_aip)
        call generate_a2(xyz, w, n, 0.1907581241803167e-3_aip)
        call generate_a3(xyz, w, n, 0.1901059546737578e-3_aip)
        call generate_bk(xyz, w, 0.1182361662400277e-1_aip, n, 0.3926424538919212e-4_aip)
        call generate_bk(xyz, w, 0.3062145009138958e-1_aip, n, 0.6667905467294382e-4_aip)
        call generate_bk(xyz, w, 0.5329794036834243e-1_aip, n, 0.8868891315019135e-4_aip)
        call generate_bk(xyz, w, 0.7848165532862220e-1_aip, n, 0.1066306000958872e-3_aip)
        call generate_bk(xyz, w, 0.1054038157636201_aip, n, 0.1214506743336128e-3_aip)
        call generate_bk(xyz, w, 0.1335577797766211_aip, n, 0.1338054681640871e-3_aip)
        call generate_bk(xyz, w, 0.1625769955502252_aip, n, 0.1441677023628504e-3_aip)
        call generate_bk(xyz, w, 0.1921787193412792_aip, n, 0.1528880200826557e-3_aip)
        call generate_bk(xyz, w, 0.2221340534690548_aip, n, 0.1602330623773609e-3_aip)
        call generate_bk(xyz, w, 0.2522504912791132_aip, n, 0.1664102653445244e-3_aip)
        call generate_bk(xyz, w, 0.2823610860679697_aip, n, 0.1715845854011323e-3_aip)
        call generate_bk(xyz, w, 0.3123173966267560_aip, n, 0.1758901000133069e-3_aip)
        call generate_bk(xyz, w, 0.3419847036953789_aip, n, 0.1794382485256736e-3_aip)
        call generate_bk(xyz, w, 0.3712386456999758_aip, n, 0.1823238106757407e-3_aip)
        call generate_bk(xyz, w, 0.3999627649876828_aip, n, 0.1846293252959976e-3_aip)
        call generate_bk(xyz, w, 0.4280466458648093_aip, n, 0.1864284079323098e-3_aip)
        call generate_bk(xyz, w, 0.4553844360185711_aip, n, 0.1877882694626914e-3_aip)
        call generate_bk(xyz, w, 0.4818736094437834_aip, n, 0.1887716321852025e-3_aip)
        call generate_bk(xyz, w, 0.5074138709260629_aip, n, 0.1894381638175673e-3_aip)
        call generate_bk(xyz, w, 0.5319061304570707_aip, n, 0.1898454899533629e-3_aip)
        call generate_bk(xyz, w, 0.5552514978677286_aip, n, 0.1900497929577815e-3_aip)
        call generate_bk(xyz, w, 0.5981009025246183_aip, n, 0.1900671501924092e-3_aip)
        call generate_bk(xyz, w, 0.6173990192228116_aip, n, 0.1899837555533510e-3_aip)
        call generate_bk(xyz, w, 0.6351365239411131_aip, n, 0.1899014113156229e-3_aip)
        call generate_bk(xyz, w, 0.6512010228227200_aip, n, 0.1898581257705106e-3_aip)
        call generate_bk(xyz, w, 0.6654758363948120_aip, n, 0.1898804756095753e-3_aip)
        call generate_bk(xyz, w, 0.6778410414853370_aip, n, 0.1899793610426402e-3_aip)
        call generate_bk(xyz, w, 0.6881760887484110_aip, n, 0.1901464554844117e-3_aip)
        call generate_bk(xyz, w, 0.6963645267094598_aip, n, 0.1903533246259542e-3_aip)
        call generate_bk(xyz, w, 0.7023010617153579_aip, n, 0.1905556158463228e-3_aip)
        call generate_bk(xyz, w, 0.7059004636628753_aip, n, 0.1907037155663528e-3_aip)
        call generate_ck(xyz, w, 0.3552470312472575e-1_aip, n, 0.5992997844249967e-4_aip)
        call generate_ck(xyz, w, 0.9151176620841283e-1_aip, n, 0.9749059382456978e-4_aip)
        call generate_ck(xyz, w, 0.1566197930068980_aip, n, 0.1241680804599158e-3_aip)
        call generate_ck(xyz, w, 0.2265467599271907_aip, n, 0.1437626154299360e-3_aip)
        call generate_ck(xyz, w, 0.2988242318581361_aip, n, 0.1584200054793902e-3_aip)
        call generate_ck(xyz, w, 0.3717482419703886_aip, n, 0.1694436550982744e-3_aip)
        call generate_ck(xyz, w, 0.4440094491758889_aip, n, 0.1776617014018108e-3_aip)
        call generate_ck(xyz, w, 0.5145337096756642_aip, n, 0.1836132434440077e-3_aip)
        call generate_ck(xyz, w, 0.5824053672860230_aip, n, 0.1876494727075983e-3_aip)
        call generate_ck(xyz, w, 0.6468283961043370_aip, n, 0.1899906535336482e-3_aip)
        call generate_dk(xyz, w, 0.6095964259104373e-1_aip, 0.1787828275342931e-1_aip, n, 0.8143252820767350e-4_aip)
        call generate_dk(xyz, w, 0.8811962270959388e-1_aip, 0.3953888740792096e-1_aip, n, 0.9998859890887728e-4_aip)
        call generate_dk(xyz, w, 0.1165936722428831_aip, 0.6378121797722990e-1_aip, n, 0.1156199403068359e-3_aip)
        call generate_dk(xyz, w, 0.1460232857031785_aip, 0.8985890813745037e-1_aip, n, 0.1287632092635513e-3_aip)
        call generate_dk(xyz, w, 0.1761197110181755_aip, 0.1172606510576162_aip, n, 0.1398378643365139e-3_aip)
        call generate_dk(xyz, w, 0.2066471190463718_aip, 0.1456102876970995_aip, n, 0.1491876468417391e-3_aip)
        call generate_dk(xyz, w, 0.2374076026328152_aip, 0.1746153823011775_aip, n, 0.1570855679175456e-3_aip)
        call generate_dk(xyz, w, 0.2682305474337051_aip, 0.2040383070295584_aip, n, 0.1637483948103775e-3_aip)
        call generate_dk(xyz, w, 0.2989653312142369_aip, 0.2336788634003698_aip, n, 0.1693500566632843e-3_aip)
        call generate_dk(xyz, w, 0.3294762752772209_aip, 0.2633632752654219_aip, n, 0.1740322769393633e-3_aip)
        call generate_dk(xyz, w, 0.3596390887276086_aip, 0.2929369098051601_aip, n, 0.1779126637278296e-3_aip)
        call generate_dk(xyz, w, 0.3893383046398812_aip, 0.3222592785275512_aip, n, 0.1810908108835412e-3_aip)
        call generate_dk(xyz, w, 0.4184653789358347_aip, 0.3512004791195743_aip, n, 0.1836529132600190e-3_aip)
        call generate_dk(xyz, w, 0.4469172319076166_aip, 0.3796385677684537_aip, n, 0.1856752841777379e-3_aip)
        call generate_dk(xyz, w, 0.4745950813276976_aip, 0.4074575378263879_aip, n, 0.1872270566606832e-3_aip)
        call generate_dk(xyz, w, 0.5014034601410262_aip, 0.4345456906027828_aip, n, 0.1883722645591307e-3_aip)
        call generate_dk(xyz, w, 0.5272493404551239_aip, 0.4607942515205134_aip, n, 0.1891714324525297e-3_aip)
        call generate_dk(xyz, w, 0.5520413051846366_aip, 0.4860961284181720_aip, n, 0.1896827480450146e-3_aip)
        call generate_dk(xyz, w, 0.5756887237503077_aip, 0.5103447395342790_aip, n, 0.1899628417059528e-3_aip)
        call generate_dk(xyz, w, 0.1225039430588352_aip, 0.2136455922655793e-1_aip, n, 0.1123301829001669e-3_aip)
        call generate_dk(xyz, w, 0.1539113217321372_aip, 0.4520926166137188e-1_aip, n, 0.1253698826711277e-3_aip)
        call generate_dk(xyz, w, 0.1856213098637712_aip, 0.7086468177864818e-1_aip, n, 0.1366266117678531e-3_aip)
        call generate_dk(xyz, w, 0.2174998728035131_aip, 0.9785239488772918e-1_aip, n, 0.1462736856106918e-3_aip)
        call generate_dk(xyz, w, 0.2494128336938330_aip, 0.1258106396267210_aip, n, 0.1545076466685412e-3_aip)
        call generate_dk(xyz, w, 0.2812321562143480_aip, 0.1544529125047001_aip, n, 0.1615096280814007e-3_aip)
        call generate_dk(xyz, w, 0.3128372276456111_aip, 0.1835433512202753_aip, n, 0.1674366639741759e-3_aip)
        call generate_dk(xyz, w, 0.3441145160177973_aip, 0.2128813258619585_aip, n, 0.1724225002437900e-3_aip)
        call generate_dk(xyz, w, 0.3749567714853510_aip, 0.2422913734880829_aip, n, 0.1765810822987288e-3_aip)
        call generate_dk(xyz, w, 0.4052621732015610_aip, 0.2716163748391453_aip, n, 0.1800104126010751e-3_aip)
        call generate_dk(xyz, w, 0.4349335453522385_aip, 0.3007127671240280_aip, n, 0.1827960437331284e-3_aip)
        call generate_dk(xyz, w, 0.4638776641524965_aip, 0.3294470677216479_aip, n, 0.1850140300716308e-3_aip)
        call generate_dk(xyz, w, 0.4920046410462687_aip, 0.3576932543699155_aip, n, 0.1867333507394938e-3_aip)
        call generate_dk(xyz, w, 0.5192273554861704_aip, 0.3853307059757764_aip, n, 0.1880178688638289e-3_aip)
        call generate_dk(xyz, w, 0.5454609081136522_aip, 0.4122425044452694_aip, n, 0.1889278925654758e-3_aip)
        call generate_dk(xyz, w, 0.5706220661424140_aip, 0.4383139587781027_aip, n, 0.1895213832507346e-3_aip)
        call generate_dk(xyz, w, 0.5946286755181518_aip, 0.4634312536300553_aip, n, 0.1898548277397420e-3_aip)
        call generate_dk(xyz, w, 0.1905370790924295_aip, 0.2371311537781979e-1_aip, n, 0.1349105935937341e-3_aip)
        call generate_dk(xyz, w, 0.2242518717748009_aip, 0.4917878059254806e-1_aip, n, 0.1444060068369326e-3_aip)
        call generate_dk(xyz, w, 0.2577190808025936_aip, 0.7595498960495142e-1_aip, n, 0.1526797390930008e-3_aip)
        call generate_dk(xyz, w, 0.2908724534927187_aip, 0.1036991083191100_aip, n, 0.1598208771406474e-3_aip)
        call generate_dk(xyz, w, 0.3236354020056219_aip, 0.1321348584450234_aip, n, 0.1659354368615331e-3_aip)
        call generate_dk(xyz, w, 0.3559267359304543_aip, 0.1610316571314789_aip, n, 0.1711279910946440e-3_aip)
        call generate_dk(xyz, w, 0.3876637123676956_aip, 0.1901912080395707_aip, n, 0.1754952725601440e-3_aip)
        call generate_dk(xyz, w, 0.4187636705218842_aip, 0.2194384950137950_aip, n, 0.1791247850802529e-3_aip)
        call generate_dk(xyz, w, 0.4491449019883107_aip, 0.2486155334763858_aip, n, 0.1820954300877716e-3_aip)
        call generate_dk(xyz, w, 0.4787270932425445_aip, 0.2775768931812335_aip, n, 0.1844788524548449e-3_aip)
        call generate_dk(xyz, w, 0.5074315153055574_aip, 0.3061863786591120_aip, n, 0.1863409481706220e-3_aip)
        call generate_dk(xyz, w, 0.5351810507738336_aip, 0.3343144718152556_aip, n, 0.1877433008795068e-3_aip)
        call generate_dk(xyz, w, 0.5619001025975381_aip, 0.3618362729028427_aip, n, 0.1887444543705232e-3_aip)
        call generate_dk(xyz, w, 0.5875144035268046_aip, 0.3886297583620408_aip, n, 0.1894009829375006e-3_aip)
        call generate_dk(xyz, w, 0.6119507308734495_aip, 0.4145742277792031_aip, n, 0.1897683345035198e-3_aip)
        call generate_dk(xyz, w, 0.2619733870119463_aip, 0.2540047186389353e-1_aip, n, 0.1517327037467653e-3_aip)
        call generate_dk(xyz, w, 0.2968149743237949_aip, 0.5208107018543989e-1_aip, n, 0.1587740557483543e-3_aip)
        call generate_dk(xyz, w, 0.3310451504860488_aip, 0.7971828470885599e-1_aip, n, 0.1649093382274097e-3_aip)
        call generate_dk(xyz, w, 0.3646215567376676_aip, 0.1080465999177927_aip, n, 0.1701915216193265e-3_aip)
        call generate_dk(xyz, w, 0.3974916785279360_aip, 0.1368413849366629_aip, n, 0.1746847753144065e-3_aip)
        call generate_dk(xyz, w, 0.4295967403772029_aip, 0.1659073184763559_aip, n, 0.1784555512007570e-3_aip)
        call generate_dk(xyz, w, 0.4608742854473447_aip, 0.1950703730454614_aip, n, 0.1815687562112174e-3_aip)
        call generate_dk(xyz, w, 0.4912598858949903_aip, 0.2241721144376724_aip, n, 0.1840864370663302e-3_aip)
        call generate_dk(xyz, w, 0.5206882758945558_aip, 0.2530655255406489_aip, n, 0.1860676785390006e-3_aip)
        call generate_dk(xyz, w, 0.5490940914019819_aip, 0.2816118409731066_aip, n, 0.1875690583743703e-3_aip)
        call generate_dk(xyz, w, 0.5764123302025542_aip, 0.3096780504593238_aip, n, 0.1886453236347225e-3_aip)
        call generate_dk(xyz, w, 0.6025786004213506_aip, 0.3371348366394987_aip, n, 0.1893501123329645e-3_aip)
        call generate_dk(xyz, w, 0.6275291964794956_aip, 0.3638547827694396_aip, n, 0.1897366184519868e-3_aip)
        call generate_dk(xyz, w, 0.3348189479861771_aip, 0.2664841935537443e-1_aip, n, 0.1643908815152736e-3_aip)
        call generate_dk(xyz, w, 0.3699515545855295_aip, 0.5424000066843495e-1_aip, n, 0.1696300350907768e-3_aip)
        call generate_dk(xyz, w, 0.4042003071474669_aip, 0.8251992715430854e-1_aip, n, 0.1741553103844483e-3_aip)
        call generate_dk(xyz, w, 0.4375320100182624_aip, 0.1112695182483710_aip, n, 0.1780015282386092e-3_aip)
        call generate_dk(xyz, w, 0.4699054490335947_aip, 0.1402964116467816_aip, n, 0.1812116787077125e-3_aip)
        call generate_dk(xyz, w, 0.5012739879431952_aip, 0.1694275117584291_aip, n, 0.1838323158085421e-3_aip)
        call generate_dk(xyz, w, 0.5315874883754966_aip, 0.1985038235312689_aip, n, 0.1859113119837737e-3_aip)
        call generate_dk(xyz, w, 0.5607937109622117_aip, 0.2273765660020893_aip, n, 0.1874969220221698e-3_aip)
        call generate_dk(xyz, w, 0.5888393223495521_aip, 0.2559041492849764_aip, n, 0.1886375612681076e-3_aip)
        call generate_dk(xyz, w, 0.6156705979160163_aip, 0.2839497251976899_aip, n, 0.1893819575809276e-3_aip)
        call generate_dk(xyz, w, 0.6412338809078123_aip, 0.3113791060500690_aip, n, 0.1897794748256767e-3_aip)
        call generate_dk(xyz, w, 0.4076051259257167_aip, 0.2757792290858463e-1_aip, n, 0.1738963926584846e-3_aip)
        call generate_dk(xyz, w, 0.4423788125791520_aip, 0.5584136834984293e-1_aip, n, 0.1777442359873466e-3_aip)
        call generate_dk(xyz, w, 0.4760480917328258_aip, 0.8457772087727143e-1_aip, n, 0.1810010815068719e-3_aip)
        call generate_dk(xyz, w, 0.5085838725946297_aip, 0.1135975846359248_aip, n, 0.1836920318248129e-3_aip)
        call generate_dk(xyz, w, 0.5399513637391218_aip, 0.1427286904765053_aip, n, 0.1858489473214328e-3_aip)
        call generate_dk(xyz, w, 0.5701118433636380_aip, 0.1718112740057635_aip, n, 0.1875079342496592e-3_aip)
        call generate_dk(xyz, w, 0.5990240530606021_aip, 0.2006944855985351_aip, n, 0.1887080239102310e-3_aip)
        call generate_dk(xyz, w, 0.6266452685139695_aip, 0.2292335090598907_aip, n, 0.1894905752176822e-3_aip)
        call generate_dk(xyz, w, 0.6529320971415942_aip, 0.2572871512353714_aip, n, 0.1898991061200695e-3_aip)
        call generate_dk(xyz, w, 0.4791583834610126_aip, 0.2826094197735932e-1_aip, n, 0.1809065016458791e-3_aip)
        call generate_dk(xyz, w, 0.5130373952796940_aip, 0.5699871359683649e-1_aip, n, 0.1836297121596799e-3_aip)
        call generate_dk(xyz, w, 0.5456252429628476_aip, 0.8602712528554394e-1_aip, n, 0.1858426916241869e-3_aip)
        call generate_dk(xyz, w, 0.5768956329682385_aip, 0.1151748137221281_aip, n, 0.1875654101134641e-3_aip)
        call generate_dk(xyz, w, 0.6068186944699046_aip, 0.1442811654136362_aip, n, 0.1888240751833503e-3_aip)
        call generate_dk(xyz, w, 0.6353622248024907_aip, 0.1731930321657680_aip, n, 0.1896497383866979e-3_aip)
        call generate_dk(xyz, w, 0.6624927035731797_aip, 0.2017619958756061_aip, n, 0.1900775530219121e-3_aip)
        call generate_dk(xyz, w, 0.5484933508028488_aip, 0.2874219755907391e-1_aip, n, 0.1858525041478814e-3_aip)
        call generate_dk(xyz, w, 0.5810207682142106_aip, 0.5778312123713695e-1_aip, n, 0.1876248690077947e-3_aip)
        call generate_dk(xyz, w, 0.6120955197181352_aip, 0.8695262371439526e-1_aip, n, 0.1889404439064607e-3_aip)
        call generate_dk(xyz, w, 0.6416944284294319_aip, 0.1160893767057166_aip, n, 0.1898168539265290e-3_aip)
        call generate_dk(xyz, w, 0.6697926391731260_aip, 0.1450378826743251_aip, n, 0.1902779940661772e-3_aip)
        call generate_dk(xyz, w, 0.6147594390585488_aip, 0.2904957622341456e-1_aip, n, 0.1890125641731815e-3_aip)
        call generate_dk(xyz, w, 0.6455390026356783_aip, 0.5823809152617197e-1_aip, n, 0.1899434637795751e-3_aip)
        call generate_dk(xyz, w, 0.6747258588365477_aip, 0.8740384899884715e-1_aip, n, 0.1904520856831751e-3_aip)
        call generate_dk(xyz, w, 0.6772135750395347_aip, 0.2919946135808105e-1_aip, n, 0.1905534498734563e-3_aip)

        ! Obtain the spherical angular coordinates
        ang = cartesian_2_spherical(xyz, 1.0_aip)

        ! We finally include de 4*pi term in the weight
        w(:) = 4.0_aip * pi * w(:)
    
    end subroutine compute_angular_mesh_lebedev_131
        
    !> Generate a set of points for points with 6 equivalences with generator (1,0,0) known as \(a_1\):
    !> @param[in,out] xyz - the mesh point
    !> @param[in,out] w - the weights
    !> @param[in,ouy] n - the starting index
    !> @param[in] set_weight - the weight of that set of points
    subroutine generate_a1(xyz, w, n, set_weight)

        real(aip),    intent(inout) :: xyz(:,:)
        real(aip),    intent(inout) :: w(:)
        integer(i32), intent(inout) :: n
        real(aip),    intent(in)    :: set_weight


        xyz(n,:)   = [ 1.0_aip , 0.0_aip , 0.0_aip]
        xyz(n+1,:) = [-1.0_aip , 0.0_aip , 0.0_aip]
        xyz(n+2,:) = [ 0.0_aip , 1.0_aip , 0.0_aip]
        xyz(n+3,:) = [ 0.0_aip ,-1.0_aip , 0.0_aip]
        xyz(n+4,:) = [ 0.0_aip , 0.0_aip , 1.0_aip]
        xyz(n+5,:) = [ 0.0_aip , 0.0_aip ,-1.0_aip]

        w(n:n+5)   =  set_weight

        n = n + 6

    end subroutine generate_a1

    !> Generate a set of points for points with 12 equivalences known as \(a_2\):
    !> @param[in,out] xyz - the mesh point
    !> @param[in,out] w - the weights
    !> @param[in] n - the starting index
    !> @param[in] set_weight - the weight of that set of points
    subroutine generate_a2(xyz, w, n, set_weight)

        real(aip),    intent(inout) :: xyz(:,:)
        real(aip),    intent(inout) :: w(:)
        integer(i32), intent(inout) :: n
        real(aip),    intent(in)    :: set_weight

        real(aip), parameter :: a = sqrt(0.5_aip)
        real(aip), parameter :: z = 0.0_aip

        xyz(n+0 ,:) = [  z,   a,   a]
        xyz(n+1 ,:) = [  z,   a,  -a]
        xyz(n+2 ,:) = [  z,  -a,   a]
        xyz(n+3 ,:) = [  z,  -a,  -a]
        xyz(n+4 ,:) = [  a,   z,   a]
        xyz(n+5 ,:) = [  a,   z,  -a]
        xyz(n+6 ,:) = [ -a,   z,   a]
        xyz(n+7 ,:) = [ -a,   z,  -a]
        xyz(n+8 ,:) = [  a,   a,   z]
        xyz(n+9 ,:) = [  a,  -a,   z]
        xyz(n+10,:) = [ -a,   a,   z]
        xyz(n+11,:) = [ -a,  -a,   z]

        w(n:n+11)   =  set_weight

        n = n + 12

    end subroutine generate_a2

    !> Generate a set of points for points with 8 equivalences known as \(a_3\):
    !> @param[in,out] xyz - the mesh point
    !> @param[in,out] w - the weights
    !> @param[in] n - the starting index
    !> @param[in] set_weight - the weight of that set of points
    subroutine generate_a3(xyz, w, n, set_weight)

        real(aip),    intent(inout) :: xyz(:,:)
        real(aip),    intent(inout) :: w(:)
        integer(i32), intent(inout) :: n
        real(aip),    intent(in)    :: set_weight

        real(aip), parameter :: a = 1.0_aip / sqrt(3.0_aip)

        xyz(n   ,:) = [  a,  a,  a]
        xyz(n+1 ,:) = [  a,  a, -a]
        xyz(n+2 ,:) = [  a, -a,  a]
        xyz(n+3 ,:) = [  a, -a, -a]
        xyz(n+4 ,:) = [ -a,  a,  a]
        xyz(n+5 ,:) = [ -a,  a, -a]
        xyz(n+6 ,:) = [ -a, -a,  a]
        xyz(n+7 ,:) = [ -a, -a, -a]

        w(n:n+7)   =  set_weight

        n = n + 8

    end subroutine generate_a3

    !> Generate a set of points for points with 24 equivalences known as \(b_k\):
    !> with points \([l_k,l_k,m_k]\) and a constraint \(2l_k^2 + m_k^2 = 1\)
    !> @param[in,out] xyz - the mesh point
    !> @param[in,out] w - the weights
    !> @param[in] lk - \(l_k\) element of the point set
    !> @param[in] n - the starting index
    !> @param[in] set_weight - the weight of that set of points
    subroutine generate_bk(xyz, w, lk, n, set_weight)

        real(aip),    intent(inout) :: xyz(:,:)
        real(aip),    intent(inout) :: w(:)
        real(aip),    intent(in)    :: lk
        integer(i32), intent(inout) :: n
        real(aip),    intent(in)    :: set_weight

        real(aip) :: mk

        mk = sqrt(1.0_aip-2.0_aip *lk**2)

        xyz(n,  :)  =  [ lk,   lk,   mk ]
        xyz(n+1,:)  =  [ lk,   lk,  -mk ]
        xyz(n+2,:)  =  [ lk,  -lk,   mk ]
        xyz(n+3,:)  =  [ lk,  -lk,  -mk ]
        xyz(n+4,:)  =  [-lk,   lk,   mk ]
        xyz(n+5,:)  =  [-lk,   lk,  -mk ]
        xyz(n+6,:)  =  [-lk,  -lk,   mk ]
        xyz(n+7,:)  =  [-lk,  -lk,  -mk ]
        xyz(n+8,:)  =  [ lk,   mk,   lk ]
        xyz(n+9,:)  =  [ lk,  -mk,   lk ]
        xyz(n+10,:) =  [ lk,   mk,  -lk ]
        xyz(n+11,:) =  [ lk,  -mk,  -lk ]
        xyz(n+12,:) =  [-lk,   mk,   lk ]
        xyz(n+13,:) =  [-lk,  -mk,   lk ]
        xyz(n+14,:) =  [-lk,   mk,  -lk ]
        xyz(n+15,:) =  [-lk,  -mk,  -lk ]
        xyz(n+16,:) =  [ mk,   lk,   lk ]
        xyz(n+17,:) =  [-mk,   lk,   lk ]
        xyz(n+18,:) =  [ mk,   lk,  -lk ]
        xyz(n+19,:) =  [-mk,   lk,  -lk ]
        xyz(n+20,:) =  [ mk,  -lk,   lk ]
        xyz(n+21,:) =  [-mk,  -lk,   lk ]
        xyz(n+22,:) =  [ mk,  -lk,  -lk ]
        xyz(n+23,:) =  [-mk,  -lk,  -lk ]

        w(n:n+23)   =  set_weight

        n = n + 24

    end subroutine generate_bk

    !> Generate a set of points for points with 24 equivalences known as \(c_k\):
    !> with points \([p_k,q_k,0]\) and a constraint \(p_k^2 + q_k^2 = 1\)
    !> @param[in,out] xyz - the mesh point
    !> @param[in,out] w - the weights
    !> @param[in] pk - \(p_k\) element of the point set
    !> @param[in] n - the starting index
    !> @param[in] set_weight - the weight of that set of points
    subroutine generate_ck(xyz, w, pk, n, set_weight)

        real(aip),    intent(inout) :: xyz(:,:)
        real(aip),    intent(inout) :: w(:)
        real(aip),    intent(in)    :: pk
        integer(i32), intent(inout) :: n
        real(aip),    intent(in)    :: set_weight

        real(aip) :: qk
        real(aip), parameter :: z = 0.0_aip

        qk = sqrt(1.0_aip-pk**2)

        xyz(n,  :)  = [   pk,   qk,    z]
        xyz(n+1,:)  = [   pk,  -qk,    z]
        xyz(n+2,:)  = [  -pk,   qk,    z]
        xyz(n+3,:)  = [  -pk,  -qk,    z]
        xyz(n+4,:)  = [   qk,   pk,    z]
        xyz(n+5,:)  = [   qk,  -pk,    z]
        xyz(n+6,:)  = [  -qk,   pk,    z]
        xyz(n+7,:)  = [  -qk,  -pk,    z]
        xyz(n+8,:)  = [   pk,    z,   qk]
        xyz(n+9,:)  = [   pk,    z,  -qk]
        xyz(n+10,:) = [  -pk,    z,   qk]
        xyz(n+11,:) = [  -pk,    z,  -qk]
        xyz(n+12,:) = [   qk,    z,   pk]
        xyz(n+13,:) = [   qk,    z,  -pk]
        xyz(n+14,:) = [  -qk,    z,   pk]
        xyz(n+15,:) = [  -qk,    z,  -pk]
        xyz(n+16,:) = [    z,   pk,   qk]
        xyz(n+17,:) = [    z,   pk,  -qk]
        xyz(n+18,:) = [    z,  -pk,   qk]
        xyz(n+19,:) = [    z,  -pk,  -qk]
        xyz(n+20,:) = [    z,   qk,   pk]
        xyz(n+21,:) = [    z,   qk,  -pk]
        xyz(n+22,:) = [    z,  -qk,   pk]
        xyz(n+23,:) = [    z,  -qk,  -pk]

        w(n:n+23)   =  set_weight

        n = n + 24

    end subroutine generate_ck

    !> Generate a set of points for points with 48 equivalences known as \(d_k\):
    !> with points \([r_k,S_k,W_k]\) and a constraint \(r_k^2 + S_k^2 + W_k^2 = 1\)
    !> @param[in,out] xyz - the mesh point
    !> @param[in,out] w - the weights
    !> @param[in] rk - \(r_k\) element of the point set
    !> @param[in] Sk - \(S_k\) element of the point set
    !> @param[in] n - the starting index
    !> @param[in] set_weight - the weight of that set of points
    subroutine generate_dk(xyz, w, rk, Sk, n, set_weight)

        real(aip),    intent(inout) :: xyz(:,:)
        real(aip),    intent(inout) :: w(:)
        real(aip),    intent(in)    :: rk
        real(aip),    intent(in)    :: Sk
        integer(i32), intent(inout) :: n
        real(aip),    intent(in)    :: set_weight

        real(aip) :: Wk

        Wk = sqrt(1.0_aip-rk**2-Sk**2)

        xyz(n,:)    = [  rk,   Sk,   Wk ]
        xyz(n+1,:)  = [  rk,   Sk,  -Wk ]
        xyz(n+2,:)  = [  rk,  -Sk,   Wk ]
        xyz(n+3,:)  = [  rk,  -Sk,  -Wk ]
        xyz(n+4,:)  = [ -rk,   Sk,   Wk ]
        xyz(n+5,:)  = [ -rk,   Sk,  -Wk ]
        xyz(n+6,:)  = [ -rk,  -Sk,   Wk ]
        xyz(n+7,:)  = [ -rk,  -Sk,  -Wk ]
        xyz(n+8,:)  = [  rk,   Wk,   Sk ]
        xyz(n+9,:)  = [  rk,   Wk,  -Sk ]
        xyz(n+10,:) = [  rk,  -Wk,   Sk ]
        xyz(n+11,:) = [  rk,  -Wk,  -Sk ]
        xyz(n+12,:) = [ -rk,   Wk,   Sk ]
        xyz(n+13,:) = [ -rk,   Wk,  -Sk ]
        xyz(n+14,:) = [ -rk,  -Wk,   Sk ]
        xyz(n+15,:) = [ -rk,  -Wk,  -Sk ]
        xyz(n+16,:) = [  Sk,   rk,   Wk ]
        xyz(n+17,:) = [  Sk,   rk,  -Wk ]
        xyz(n+18,:) = [  Sk,  -rk,   Wk ]
        xyz(n+19,:) = [  Sk,  -rk,  -Wk ]
        xyz(n+20,:) = [ -Sk,   rk,   Wk ]
        xyz(n+21,:) = [ -Sk,   rk,  -Wk ]
        xyz(n+22,:) = [ -Sk,  -rk,   Wk ]
        xyz(n+23,:) = [ -Sk,  -rk,  -Wk ]
        xyz(n+24,:) = [  Sk,   Wk,   rk ]
        xyz(n+25,:) = [  Sk,   Wk,  -rk ]
        xyz(n+26,:) = [  Sk,  -Wk,   rk ]
        xyz(n+27,:) = [  Sk,  -Wk,  -rk ]
        xyz(n+28,:) = [ -Sk,   Wk,   rk ]
        xyz(n+29,:) = [ -Sk,   Wk,  -rk ]
        xyz(n+30,:) = [ -Sk,  -Wk,   rk ]
        xyz(n+31,:) = [ -Sk,  -Wk,  -rk ]
        xyz(n+32,:) = [  Wk,   rk,   Sk ]
        xyz(n+33,:) = [  Wk,   rk,  -Sk ]
        xyz(n+34,:) = [  Wk,  -rk,   Sk ]
        xyz(n+35,:) = [  Wk,  -rk,  -Sk ]
        xyz(n+36,:) = [ -Wk,   rk,   Sk ]
        xyz(n+37,:) = [ -Wk,   rk,  -Sk ]
        xyz(n+38,:) = [ -Wk,  -rk,   Sk ]
        xyz(n+39,:) = [ -Wk,  -rk,  -Sk ]
        xyz(n+40,:) = [  Wk,   Sk,   rk ]
        xyz(n+41,:) = [  Wk,   Sk,  -rk ]
        xyz(n+42,:) = [  Wk,  -Sk,   rk ]
        xyz(n+43,:) = [  Wk,  -Sk,  -rk ]
        xyz(n+44,:) = [ -Wk,   Sk,   rk ]
        xyz(n+45,:) = [ -Wk,   Sk,  -rk ]
        xyz(n+46,:) = [ -Wk,  -Sk,   rk ]
        xyz(n+47,:) = [ -Wk,  -Sk,  -rk ]

        w(n:n+47)   =  set_weight

        n = n + 48

    end subroutine generate_dk

end module idiel_sph_quadrature
