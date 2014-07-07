Subroutine construct_t2 (irep, sym_rt2)
!
use mod_symmetry
use raman_symmetry
!
Implicit None
! arguments
integer, intent(in) :: irep
complex(8), intent(out) :: sym_rt2(3, 3, 3, 3)
! local variables
Integer :: isym, i, j, iop1, iop2, jop1, jop2
Real(8) :: sc(3, 3)
complex(8) :: s(3, 3)
character(4) :: dm(3,3), dmt(3,3)
real(8) :: t1, t1r, t1i, t2, t2r, t2i
logical :: done(3,3)
!
Do iop1 = 1, 3
   Do iop2 = 1, 3
      s (:, :) = zzero
      Do isym = 1, numsop
         sc(:, :) = sopmatc(:, :, isym)
         Do i = 1, 3
            Do j = 1, 3
               s(i, j) = s(i, j) + charact(class(isym), irep) * sc(i, iop1) * sc(j, iop2)
            End Do
         End Do
      End Do
      sym_rt2(iop1, iop2, :, :) = s(:, :) / dble(numsop)
   End Do
End Do
!
done(:, :) = .false.
dm(1, :) = (/ 'a_11', 'a_12', 'a_13' /)
dm(2, :) = (/ 'a_21', 'a_22', 'a_23' /)
dm(3, :) = (/ 'a_31', 'a_32', 'a_33' /)
dmt(:, :) = dm(:, :)
!
! analyze the symmetrization matrices
do iop1 = 1, 3
  do iop2 = 1, 3
    t2r = sum(abs(dble(sym_rt2(iop1, iop2, :, :))))
    t2i = sum(abs(aimag(sym_rt2(iop1, iop2, :, :))))
    t2 = t2r**2 + t2i**2
    do jop1 = 1, 3
      do jop2 = 1, 3
        if (.not.done(iop1, iop2)) then
          t1r = sum(abs(dble(sym_rt2(iop1, iop2, :, :)) - dble(sym_rt2(jop1, jop2, :, :))))
          t1i = sum(abs(aimag(sym_rt2(iop1, iop2, :, :)) - aimag(sym_rt2(jop1, jop2, :, :))))
          t1 = t1r**2 + t1i**2
          if ((t1.lt.eps) .and. (t2.gt.eps) .and. ((iop1.ne.jop1) .or. (iop2.ne.jop2))) then
            dmt(jop1, jop2) = dm(iop1, iop2)
          end if
          done(jop1, jop2) = .true.
       end if
      end do
    end do
    if (sum(abs(sym_rt2(iop1, iop2, :, :))) .lt. eps) then
      dmt(iop1, iop2) = ' 0  '
    end if
  end do
end do
  ! output Raman tensor symmetry
open( unit=13, file='RAMAN_SYM.OUT', status='old', action='write', position='append')
Write (13, '(/," General structure of the Raman tensor in this symmetry: ",/)')
Do iop1 = 1, 3
  write(13,'("  ( ",a,"  ",a,"  ",a," )")') dmt(iop1,:)
End Do
write(13, *)
close(13)
!
End Subroutine construct_t2
