Subroutine construct_t2 (irep, lt2)
!
use mod_symmetry
use raman_symmetry
!
Implicit None
! arguments
integer, intent(in) :: irep
logical, intent(out) :: lt2(3, 3)
! local variables
Integer :: isym, i, j, iop1, iop2, jop1, jop2
real(8) :: symt2(3, 3, 3, 3)
Real(8) :: s(3, 3), sc(3, 3)
character(256) :: strvar
character(4) :: dm(3,3), dmt(3,3)
real(8) :: t1, t2
logical :: red(3,3), done(3,3)
!
Do iop1 = 1, 3
   Do iop2 = 1, 3
      s (:, :) = 0.d0
      Do isym = 1, numsop
         sc(:, :) = dble(sopmatc(:, :, isym))
         Do i = 1, 3
            Do j = 1, 3
               s(i, j) = s(i, j) + dble(charact(class(isym), irep)) * sc(i, iop1) * sc(j, iop2)
            End Do
         End Do
      End Do
      symt2(iop1, iop2, :, :) = s(:, :) / dble(numsop)
   End Do
End Do
!
red(:, :) = .false.
done(:, :) = .false.
lt2 = .true.
dm(1, :) = (/ 'a_11', 'a_12', 'a_13' /)
dm(2, :) = (/ 'a_21', 'a_22', 'a_23' /)
dm(3, :) = (/ 'a_31', 'a_32', 'a_33' /)
dmt(:, :) = dm(:, :)
!
! analyze the symmetrization matrices
do iop1 = 1, 3
  do iop2 = 1, 3
    t2 = sum(abs(symt2(iop1, iop2, :, :)))
    do jop1 = 1, 3
      do jop2 = 1, 3
        if (.not.done(iop1, iop2)) then
          t1 = sum(abs(symt2(iop1, iop2, :, :) - symt2(jop1, jop2, :, :)))
          if ((t1.lt.eps) .and. (t2.gt.eps) .and. ((iop1.ne.jop1) .or. (iop2.ne.jop2))) then
            dmt(jop1, jop2) = dm(iop1, iop2)
            red(jop1, jop2) = .true.
          end if
          done(jop1, jop2) = .true.
       end if
      end do
    end do
    if (sum(abs(symt2(iop1, iop2, :, :))) .lt. eps) then
      dmt(iop1, iop2) = ' 0  '
      red(iop1, iop2) = .true.
      lt2(iop1, iop2) = .false.
    end if
  end do
end do
  ! output the symmetrization matrices
open( unit=13, file='RAMAN_SYM.OUT', status='old', action='write', position='append')
Write (13, '(/," Structure of the Raman tensor: ",/)')
Do iop1 = 1, 3
  write(13,'("  ( ",a,"  ",a,"  ",a," )")') dmt(iop1,:)
End Do
write(13, *)
close(13)
!
End Subroutine construct_t2
