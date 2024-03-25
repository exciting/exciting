module mod_bravais_lattice
    use precision, only: wp
    use mod_kpoint, only: ivk, ivknr
    use mod_lattice, only: avec
    use modgw, only: kset, kqset
    
    implicit none
    private
    public :: bravais_lattice

    contains

    subroutine bravais_lattice(ik, bigR)
        !> k index
        integer, intent(in) :: ik
        !> Bravais lattice (R): R = k1*a1 + k2*a2 + k3*a3
        real(wp), intent(out) :: bigR(:)

        !write(31,*) ik, kqset%vkl(:,ik)
        !write(32,*) ik, ivk(:,ik)
        !write(33, '(7i3)') ik, ivknr(:,ik), kset%ngridk
        !write(34,*) avec(1,:)
        !write(34,*) avec(2,:)
        !write(34,*) avec(3,:)
        !! TODO (Manoar) check if we need to devide by kset%ngridk(:)
        bigR(:) = ivk(1,ik)*avec(:,1) + ivk(2,ik)*avec(:,2) + ivk(3,ik)*avec(:,3)
        !write(35, '(4i3, 3f16.11)') ik, ivk(:,ik), bigR
        ! bigR(1) = ivk(1,ik)*avec(1,1) + ivk(2,ik)*avec(1,2) + ivk(3,ik)*avec(1,3)
        ! bigR(2) = ivk(1,ik)*avec(2,1) + ivk(2,ik)*avec(2,2) + ivk(3,ik)*avec(2,3)
        ! bigR(3) = ivk(1,ik)*avec(3,1) + ivk(2,ik)*avec(3,2) + ivk(3,ik)*avec(3,3)
        ! write(36, '(i3, 3f16.11)') ik, bigR
    end subroutine bravais_lattice
end module mod_bravais_lattice
