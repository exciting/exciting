!> This module contains related stuff for computing, retrieving and storing Gaunt coefficients in a compact way, as
!> (mrm): Modified to make getgauntcoef callable from GPUs
module mod_gaunt_coefficients

    use precision, only: i32, dp
#include "offload.fpp"
    private
    public  :: calcgauntcoef, getgauntcoef, delete_gaunt_coefficients, epsangint

    real(dp), parameter :: epsangint = 1.0e-8_dp

    ! Gaunt's coefficients
    real(dp), allocatable :: gauntcoef(:)

    ! Declare here all elements for device exposure
    OMP_OFFLOAD declare target(gauntcoef)

contains
    
!---------------------------------------------------------------------------
    subroutine delete_gaunt_coefficients
        if (allocated(gauntcoef)) then 
          OMP_OFFLOAD target exit data map(delete: gauntcoef)
          deallocate(gauntcoef)
        end if
    end subroutine
    
!-------------------------------------------------------------------------------
!BOP
!!ROUTINE: calcgauntcoef
!!INTERFACE:
!
    subroutine calcgauntcoef(maxj)
      use wigner3j_symbol, only: gaunt_yyy
!
!!DESCRIPTION:
!      
! This subroutine calculates the gaunt coefficients:
! 
! \begin{equation}
! \mathcal{G}^{LM}_{ll',mm'}=\int\limits_0^{2\pi}{\int\limits_0^{\pi}{%
! Y_{lm}(\theta,\phi) Y_{l'm'}(\theta,\phi)
! Y^*_{LM}(\theta,\phi)\sin(\theta)d\theta}d\phi}
! \end{equation}
! 
! for $l$ and $l' = 0,1,$... \verb"maxj". The integral is done numerically
! on a special grid (see wigner3j\_symbol.f90).
! The values are calculated only for $l\ge l'$ and $m' \ge 0$. 
! 
! The storage is optimized by saving the values in a vector
! (\verb"gauntcoef") only for those coefficients that are different from zero.
! The size of the vector is:
!\begin{equation}
!n=\tfrac{1}{60}(l_{max}+1)(l_{max}+2)(l_{max}+3)(16l_{max}^2+29l_{max}+10)
!\end{equation}
!
!The gaunt coefficient $\mathcal{G}^{LM}_{ll',mm'}$ can be accesed directly at
!\verb"gauntcoef(i)" by applying the function:
!
!\begin{equation}
!\begin{aligned}
!i=&\tfrac{1}{60}(16l^2-3l-3)(l+2)(l+1)l+\tfrac{1}{3}ll'(l'+1)(4l'-1)+%
!\tfrac{1}{6}l'(l'-1)(4l'+7)+\\
!&(2l+1)(l'+1)(L-l+l')+(l'+1)(m+l)+m'+l'+1 
!\end{aligned}
!\end{equation}
!
!of course, $M=m+m'$ is already taken into account. For the cases $l'>l$
!and $m'<0$ see the \verb"getcgcoef" subroutine.
!
!!INPUT PARAMETERS:
        implicit none
        integer(i32), intent(in) :: maxj

!!LOCAL VARIABLES:
        integer(i32) :: i
        integer(i32) :: l1, l2,  l3
        integer(i32) :: m1, m2, m3
        integer(i32) :: ntot

!!REVISION HISTORY:
!
! Created: Apr. 2004 by RGA
! Last modified: May 21st. 2004 by RGA
! Adapted: Nov 2013 by DIN
!
!EOP
!BOC
        ntot = (maxj+1)*(maxj+2)*(maxj+3)*(16*maxj*maxj+29*maxj+10)/60
        if (allocated(gauntcoef)) then
          OMP_OFFLOAD target exit data map(delete: gauntcoef)
          deallocate(gauntcoef)
        end if
        allocate(gauntcoef(ntot), source=0.0_dp)
        i = 0
        do l1 = 0, maxj
          do l2 = 0, l1
            do l3 = l1-l2, l1+l2
              do m1 = -l1, l1
                do m2 = 0, l2
                  m3 = m1+m2
                  i = i+1
                  if (mod(l1+l2+l3,2)==0 .and. abs(m3)<=l3) then
                      gauntcoef(i) = gaunt_yyy(l3,l1,l2,m3,m1,m2)
                  endif
                end do
              end do
            end do
          end do
        end do

        OMP_OFFLOAD target enter data map(always, to: gauntcoef)

        return
    end subroutine
!EOC

!-------------------------------------------------------------------------------    
!BOP
!!ROUTINE: getgauntcoef
!!INTERFACE:
!
    real(dp) function getgauntcoef(l1,l2,l3,m1,m2)
!
!!DESCRIPTION:
!
!This function gets the gaunt coefficient $G^{l3,m1+m2}_{l1,l2,m1,m2}$
!from the vector \verb"gauntcoef" by:
!\begin{equation}
!G^{l3,m1+m2}_{l1,l2,m1,m2}=\alpha \verb"cgcoef"(i)
!\end{equation}
!calculating i by:
!\begin{equation}
!\begin{aligned}
!i=&\tfrac{1}{60}(16l^2-3l-3)(l+2)(l+1)l+\tfrac{1}{3}ll'(l'+1)(4l'-1)+%
!\tfrac{1}{6}l'(l'-1)(4l'+7)+\\
!&+(2l+1)(l'+1)(L-l+l')+(l'+1)(m+l)+m'+l'+1 
!\end{aligned}
!\end{equation}
!where
!
!\begin{subequations}
!\begin{align}
!l&=l1 & l'&=l2 & m&=m1 & m'&=m2 & \alpha&=1 &\text{ if $l1 \ge l2$} \\
!l&=l2 & l'&=l1 & m&=m2 & m'&=m1 & \alpha&=1 &\text{ if $l1 < l2$} \\
! &    &   &    & m&=m1 & m'&=m2 & \alpha&=1 &\text{ if $m2 \ge 0$} \\
! &    &   &    & m&=-m1 & m'&=-m2 & \alpha&=(-1)^{l+l'-L} &\text{ if $m2 < 0$}
!\end{align}
!\end{subequations}
!
!!INPUT PARAMETERS:
        implicit none
        integer(i32), intent(in) :: l1
        integer(i32), intent(in) :: l2
        integer(i32), intent(in) :: l3
        integer(i32), intent(in) :: m1
        integer(i32), intent(in) :: m2

!!LOCAL VARIABLES:
        integer(i32) :: j1, j2, mj1, mj2
        integer(i32) :: par, ing
        integer(i32) :: ind1, ind2, ind3, ind4
        real(dp) :: fact
        logical :: trcond

!! For GPU aware compilation this function is compiled also for the device
        OMP_OFFLOAD declare target

!!REVISION HISTORY:
!
! Created: Apr. 2004 by RGA
! Last modified  May 21st. 2004 by RGA
! Adapted: Nov 2013 by DIN
!
!EOP
!BOC
        par = mod(abs(l1+l2-l3),2)
        fact = 1.0_dp
        trcond = (abs(m1+m2)<=l3).and.(abs(l1-l2)<=l3).and.(l1+l2>=l3)
        if (trcond) then
          if (l1<l2) then
            j1 = l2
            mj1 = m2
            j2 = l1
            mj2 = m1
          else
            j1 = l1
            mj1 = m1
            j2 = l2
            mj2 = m2
          end if
          if (mj2<0) then
            mj2 = -mj2
            mj1 = -mj1
            fact = (-2.0_dp*par + 1.0_dp)
          end if
          ind1 = (16*j1*j1-3*j1-3)*(j1+2)*(j1+1)*j1/60
          ind2 = j1*j2*(j2+1)*(4*j2-1)/3
          ind3 = j2*(j2-1)*(4*j2+7)/6
          ind4 = (2*j1+1)*(j2+1)*(l3-j1+j2)
          ing = ind1+ind2+ind3+ind4+(j2+1)*(mj1+j1)+mj2+j2+1
          getgauntcoef = fact*gauntcoef(ing)
        else
          getgauntcoef = 0.0_dp
        endif
        return
    end function
!EOC

end module
