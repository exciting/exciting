!
!-----------------------------------------------------------------------------80.
! Copyright (C) 2013- exciting team
! This file is distributed under the terms of the GNU General Public License.
! Created       on 15-10-2013 Pasquale Pavone (exciting team)
! Modified      on 15-11-2013 Pasquale Pavone (exciting team)
! Last modified on 14-06-2014 Pasquale Pavone (exciting team)
!-----------------------------------------------------------------------------80
!
subroutine writetorque(fnum)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
! local variables
    integer :: is, ia, ias, i
    real(8) :: torque(3), cmpos(3), totmass

    totmass = 0.0
    cmpos(:) = 0.0
    torque(:) = 0.0
    atposcm(:,:,:) = 0.0

    call cluster

    do is = 1, nspecies
        do ia = 1, natoms (is)
            totmass = totmass + spmass(is)
            cmpos(:) = cmpos(:) + spmass(is)*atposcm(:,ia,is)
        end do
    end do

    cmpos(:) = cmpos(:)/totmass

    write(fnum,*)
    write(fnum,'(" Center of mass ",T18": ",3F14.8,"  (cartesian)")') cmpos(:) 

    do is = 1, nspecies
        do ia = 1, natoms (is)
            ias = idxas (ia, is)
            torque(1) = torque(1) + (atposcm(2,ia,is)-cmpos(2))*forcetot(3,ias) - &
           &                        (atposcm(3,ia,is)-cmpos(3))*forcetot(2,ias)
            torque(2) = torque(2) + (atposcm(3,ia,is)-cmpos(3))*forcetot(1,ias) - &
           &                        (atposcm(1,ia,is)-cmpos(1))*forcetot(3,ias)
            torque(3) = torque(3) + (atposcm(1,ia,is)-cmpos(1))*forcetot(2,ias) - &
           &                        (atposcm(2,ia,is)-cmpos(2))*forcetot(1,ias)
        end do
    end do

    write(fnum,'(" Total torque ",T18": ",3F14.8,"  (cartesian)")') torque(:) 
!
contains
!    
!-----------------------------------------------------------------------------80
!
    subroutine cluster
        use modinput
        use modmain

        implicit none
        integer :: i1, i2, i3
        real(8) :: d, dmin, v(3)
        real(8) :: r3dist

        external r3dist

        do is = 1, nspecies
            do ia = 1, natoms(is)
                dmin = 1.d8
                do i1 = -1, 1
                    do i2 = -1, 1
                        do i3 = -1, 1
                            v(:) = dble(i1) * input%structure%crystal%basevect(:,1) + &
                           &       dble(i2) * input%structure%crystal%basevect(:,2) + &
                           &       dble(i3) * input%structure%crystal%basevect(:,3) + &
                           &       atposc(:,ia,is)
                            d = r3dist(atposc(:,1,1), v)
                            if (d .lt. dmin) then
                                dmin = d
                                atposcm(:,ia,is) = v(:)
                            end if
                        end do
                    end do
                end do
            end do
        end do

        return
    end subroutine cluster
    
end subroutine writetorque

