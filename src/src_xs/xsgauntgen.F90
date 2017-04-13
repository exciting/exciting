! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_xsgauntgen

  implicit none

contains

  subroutine xsgauntgen(lmax1, lmax2, lmax3)
    use mod_muffin_tin, only: idxlm
    use modxs, only: xsgnt
    implicit none

    ! Arguments
    integer, intent(in) :: lmax1, lmax2, lmax3

    ! Local variables
    integer :: l1, l2, l3, m1, m2, m3
    integer :: lm1, lm2, lm3, lmmax1, lmmax2, lmmax3

    ! External functions
    real(8), external :: gaunt

    ! Allocate and generate complex gaunt coefficient array
    lmmax1 = (lmax1+1) ** 2
    lmmax2 = (lmax2+1) ** 2
    lmmax3 = (lmax3+1) ** 2

    if(allocated(xsgnt)) deallocate(xsgnt)
    allocate(xsgnt(lmmax1, lmmax2, lmmax3))

    do l1 = 0, lmax1
      do m1 = - l1, l1
        lm1 = idxlm(l1, m1)
        do l2 = 0, lmax2
          do m2 = - l2, l2
            lm2 = idxlm(l2, m2)
            do l3 = 0, lmax3
              do m3 = - l3, l3
                lm3 = idxlm(l3, m3)
                xsgnt(lm1, lm2, lm3) = gaunt(l1, l2, l3, m1, m2, m3)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine xsgauntgen

  subroutine xasgauntgen (lmax2, lmax3)
    use modxs, only: xsgntou, xsgntuo, xsgntoo
    use modxas, only: nxas, lxas, xasstart, xasstop, preml, mj2ml
    use mod_muffin_tin, only: idxlm
    Implicit None
    ! arguments
    Integer, Intent (In) :: lmax2, lmax3
    ! local variables
    Integer :: n1, n2, l2, l3, m2, m3, lm1, lm2, lm3, &
      & lmmax2, lmmax3 
    Real (8) :: gaunt, prefac
    External :: gaunt
    ! allocate and generate complex Gaunt coefficient array
    lmmax2 = (lmax2+1) ** 2
    lmmax3 = (lmax3+1) ** 2
    prefac=1.0d0/sqrt(2.0d0)
    If (allocated(xsgntou)) deallocate (xsgntou)
    If (allocated(xsgntuo)) deallocate (xsgntuo)
    If  (allocated(xsgntoo)) deallocate (xsgntoo)
    Allocate (xsgntou(nxas,lmmax2, lmmax3))
    Allocate (xsgntuo(nxas,lmmax2, lmmax3))
    Allocate (xsgntoo(nxas,lmmax2, nxas))
    Do n1 = 1, nxas             
      Do l2 = 0, lmax2
        Do m2 = - l2, l2
          lm2 = idxlm (l2, m2)
          Do l3 = 0, lmax3
            Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)
              !resonant part
              xsgntou (n1, lm2, lm3) = prefac*preml(n1+xasstart-1,1)*gaunt (lxas, l2, l3, &
                & mj2ml(n1+xasstart-1,1), m2, m3)+prefac*preml(n1+xasstart-1,2)*gaunt(lxas, l2, l3,     &
                & mj2ml(n1+xasstart-1,2), m2, m3)
              ! anti-resonant
              xsgntuo (n1, lm2, lm3) = prefac*preml(n1+xasstart-1,1)*gaunt (l3, l2, lxas, &
                & m3, m2, mj2ml(n1+xasstart-1,1))+prefac*preml(n1+xasstart-1,2)*gaunt(l3, l2, lxas,     &
                & m3, m2, mj2ml(n1+xasstart-1,2))
            End Do
          End Do
          Do n2 = 1, nxas
            xsgntoo(n1,lm2,n2)=preml(n1+xasstart-1,1)*preml(n2+xasstart-1,1) &
              &*gaunt (lxas, l2, lxas, mj2ml(n1+xasstart-1,1), m2, mj2ml &
              &(n2+xasstart-1,1))+preml(n1+xasstart-1,2)*preml(n2+xasstart-1,2)&
              &*gaunt(lxas, l2, lxas, mj2ml(n1+xasstart-1,2), m2, mj2ml(n2+xasstart-1,2))
           End Do
        End Do
      End Do
    End Do
    end subroutine xasgauntgen
end module m_xsgauntgen
