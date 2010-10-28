
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writesymt2
! !INTERFACE:
Subroutine writesymt2
! !USES:
      Use modmain
      Use modxs
! !DESCRIPTION:
!   Outputs the symmetrization matrices for the tensor components of a rank-2
!   tensor. The tensor(-field) $t_{ij}$ in reciprocal space must be invariant under
!   coordinate transforms of the system wrt. the rotational part of the crystal
!   symmetries, so we can average:
!   $$ t_{ij}^{\rm sym} = \frac{1}{N_{\alpha}}\sum_{\alpha} \sum_{k,l}
!     \alpha_{ik}\alpha_{jl}t_{kl}. $$
!   The symmetrized tensor $t_{ij}^{\rm sym}$ can then be written as
!   $$ t_{ij}^{\rm sym} = \sum_{k,l} T_{ij,kl} t_{kl}, $$
!   with the symmetrization tensor
!   $$ T_{ij,kl} = \frac{1}{N_{\alpha}}\sum_{\alpha}\alpha_{ik}\alpha_{jl} $$
!   where $N_{\alpha}$ is the number of symmetry operations in the space group.
!   For each component $ij$ the symmetrization tensor $T_{ij,kl}$ is written as
!   a matrix in the components $kl$ to the file {\tt SYMT2.OUT}.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      real(8), parameter :: eps=1.d-8
      character(256) :: strvar
      character(4) :: dm(3,3), dmt(3,3)
      Integer :: iop1, iop2, jop1, jop2, i
      real(8) :: t1, t2
      logical :: red(3,3), done(3,3)
      red(:,:)=.false.
      done(:,:)=.false.
      dm(1,:)=(/'e_11','e_12','e_13'/)
      dm(2,:)=(/'e_21','e_22','e_23'/)
      dm(3,:)=(/'e_31','e_32','e_33'/)
      dmt(:,:)=dm(:,:)
  ! analyse the symmetrization matrices
      do iop1=1,3
        do iop2=1,3
          t2=sum(abs(symt2(iop1,iop2,:,:)))
          do jop1=1,3
            do jop2=1,3
              if (.not.done(iop1,iop2)) then
                t1=sum(abs(symt2(iop1,iop2,:,:)-symt2(jop1,jop2,:,:)))
                if ((t1.lt.eps).and.(t2.gt.eps).and.((iop1.ne.jop1).or.(iop2.ne.jop2))) then
                  dmt(jop1,jop2)=dm(iop1,iop2)
                  red(jop1,jop2)=.true.
                end if
                done(jop1,jop2)=.true.
             end if
            end do
          end do
          if (sum(abs(symt2(iop1,iop2,:,:))).lt.eps) then
            dmt(iop1,iop2)=' 0  '
            red(iop1,iop2)=.true.
          end if
        end do
      end do
  ! output the symmetrization matrices
      Open (50, File='SYMT2'//trim(filext), Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("(structure of a rank 2 tensor, e.g. the dielectric tensor, &
        &from symmetry considerations")')
      Write (50, '(" with respect to Cartesian coordinates)")')
      Write (50,*)
      Write (50, '(" Upper limit for number of independent components : ",i6)') 9-count(red)
      Write (50,*)
      Do iop1 = 1, 3
        write(50,'(" ( ",a,"  ",a,"  ",a," )")') dmt(iop1,:)
      End Do
      Write (50,*)
      Write (50, '("(symmetrization matrices are in Cartesian coordinates)")')
      Write (50,*)
      Do iop1 = 1, 3
         Do iop2 = 1, 3
            t1=sum(abs(symt2(iop1,iop2,:,:)))
            strvar=""
            if (t1 .lt. eps) strvar='," ; zero contribution"'
            Write (50, '("(", i1, ", ", i2, ")-component"'//trim(strvar)//')') iop1, iop2
            Write (50, '(3f12.8)') (symt2(iop1, iop2, i, :), i=1, 3)
            Write (50,*)
         End Do
      End Do
      Close (50)
End Subroutine writesymt2
