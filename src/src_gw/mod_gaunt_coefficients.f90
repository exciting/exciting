
module mod_gaunt_coefficients

    real(8), parameter :: epsangint=1.d-8
    
    ! Gaunt's coefficients
    real(8), allocatable :: gauntcoef(:)

contains
    
!---------------------------------------------------------------------------
    subroutine delete_gaunt_coefficients
        if (allocated(gauntcoef)) deallocate(gauntcoef)
    end subroutine
    
!-------------------------------------------------------------------------------
!BOP
!!ROUTINE: calcgauntcoef
!!INTERFACE:
!
    subroutine calcgauntcoef(maxj)
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
! on a special grid (see gaunt.f90).
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
        integer(4), intent(in) :: maxj
      
!!LOCAL VARIABLES:
        integer(4) :: i
        integer(4) :: l1, l2,  l3
        integer(4) :: m1, m2, m3 
        integer(4) :: ntot, ngrid
      
!!EXTERNAL ROUTINES: 
        real(8), external :: gaunt

!!REVISION HISTORY:
!
! Created: Apr. 2004 by RGA
! Last modified: May 21st. 2004 by RGA
! Adapted: Nov 2013 by DIN
!
!EOP
!BOC
        ntot = (maxj+1)*(maxj+2)*(maxj+3)*(16*maxj*maxj+29*maxj+10)/60
        ngrid = (4*maxj+1)*(4*maxj+1)/3
        if (allocated(gauntcoef)) deallocate(gauntcoef)
        allocate(gauntcoef(ntot))
        i = 0      
        do l1 = 0, maxj
          do l2 = 0, l1
            do l3 = l1-l2, l1+l2
              do m1 = -l1, l1
                do m2 = 0, l2
                  m3 = m1+m2
                  i = i+1
                  if (mod(l1+l2+l3,2)==0) then
                    if (iabs(m3)<=l3) then
                      gauntcoef(i) = gaunt(l3,l1,l2,m3,m1,m2)
                    else
                      gauntcoef(i) = 0.0d0
                    endif
                  else
                    gauntcoef(i) = 0.0d0
                  endif    
                end do
              end do 
            end do
          end do
        end do
        return
    end subroutine
!EOC

!-------------------------------------------------------------------------------    
!BOP
!!ROUTINE: getgauntcoef
!!INTERFACE:
!
    real(8) function getgauntcoef(l1,l2,l3,m1,m2)
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
        integer(4), intent(in) :: l1
        integer(4), intent(in) :: l2
        integer(4), intent(in) :: l3
        integer(4), intent(in) :: m1
        integer(4), intent(in) :: m2
      
!!LOCAL VARIABLES:
        integer(4) :: j1, j2, mj1, mj2
        integer(4) :: par, ing
        integer(4) :: ind1, ind2, ind3, ind4
        real(8) :: fact
        logical :: trcond

!!REVISION HISTORY:
!
! Created: Apr. 2004 by RGA
! Last modified  May 21st. 2004 by RGA
! Adapted: Nov 2013 by DIN
!
!EOP
!BOC    
        par = mod(abs(l1+l2-l3),2)    
        fact = 1.0d0
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
            fact = (-2.0d0*par+1.0d0)
          end if
          ind1 = (16*j1*j1-3*j1-3)*(j1+2)*(j1+1)*j1/60
          ind2 = j1*j2*(j2+1)*(4*j2-1)/3
          ind3 = j2*(j2-1)*(4*j2+7)/6
          ind4 = (2*j1+1)*(j2+1)*(l3-j1+j2)
          ing = ind1+ind2+ind3+ind4+(j2+1)*(mj1+j1)+mj2+j2+1 
          getgauntcoef = fact*gauntcoef(ing)   
        else
          getgauntcoef = 0.0d0
        endif
        return
    end function
!EOC      

end module
