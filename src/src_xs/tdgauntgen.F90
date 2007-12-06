
module m_tdgauntgen
  implicit none
contains

  subroutine tdgauntgen(lmax1,lmax2,lmax3)
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: lmax1,lmax2,lmax3
    ! local variables
    integer :: l1,l2,l3,m1,m2,m3,lm1,lm2,lm3,lmmax1,lmmax2,lmmax3
    real(8) :: gaunt
    external :: gaunt

    ! allocate and generate complex Gaunt coefficient array
    lmmax1 = (lmax1+1)**2
    lmmax2 = (lmax2+1)**2
    lmmax3 = (lmax3+1)**2
    if (allocated(tdgnt)) deallocate(tdgnt)
    allocate(tdgnt(lmmax1,lmmax2,lmmax3))

    do l1=0,lmax1
       do m1=-l1,l1
          lm1=idxxlm(l1,m1)
          do l2=0,lmax2
             do m2=-l2,l2
                lm2=idxxlm(l2,m2)
                do l3=0,lmax3
                   do m3=-l3,l3
                      lm3=idxxlm(l3,m3)
                      tdgnt(lm1,lm2,lm3)=gaunt(l1,l2,l3,m1,m2,m3)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine tdgauntgen

end module m_tdgauntgen
