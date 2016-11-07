! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_getpemat
  use modmpi

  implicit none

contains

  !BOP
  ! !ROUTINE: getpemat
  ! !INTERFACE: getpmat
  subroutine getpemat(iq, ik, pfilnam, efilnam, m12, m34, p12, p34)
  ! !USES:
    use modinput, only: input
    use mod_constants, only: fourpi, zzero
    use modxs, only: ngq, vkl0,&
                   & istl1, istl2, istl3, istl4,&
                   & istu1, istu2, istu3, istu4,&
                   & nst1, nst2, nst3, nst4,&
                   & deou, docc12, deuo, docc21,&
                   & tscreen, xiou, xiuo, sptclg
#ifdef TETRA         
    use modtetra
#endif
    use m_getpmat
    use m_getemat
  ! !DESCRIPTION:
  ! For a given k and q point this routine reads the momentum matrix (Gamma point)
  ! and the plane wave matrix elements for the band combinations 12 and optionally
  ! 21 from file. The momentum matrix elements get renormalized by the KS transiton
  ! energies. Both momentum matrix elements and plane wave matrix elements will be 
  ! multiplied with the square root of the coulomb potential.
  !
  ! !REVISION HISTORY:
  !   Added to documentations schema. (Aurich)
  !EOP
  !BOC


    implicit none

    ! Arguments
    integer, intent(in) :: iq, ik
    character(*), intent(in) :: pfilnam, efilnam
    complex(8), optional, intent(out) :: m12(:, :, :), p34(:, :, :)
    complex(8), optional, intent(out) :: p12(:, :, :), m34(:, :, :)

    ! Local variables
    character(*), parameter :: thisnam = 'getpemat'
    real(8) :: fourpisqt
    integer :: n, igq, j, i1, i2, i3, i4
    logical :: tq0

    ! External functions
    logical, external :: tqgamma

    ! Test if iq is Gamma point
    tq0 = tqgamma(iq)

    n = ngq(iq)

    fourpisqt = sqrt(fourpi)

    ! Sanity check: For gamma point plane wave elements need
    ! to be represented by momentum matrix elements.
    if(tq0 .and. (.not. present(p12))) then
      write(*,*)
      write(*, '("Error(", a, "):&
        & Gamma q-point but momentum matrix elements not requested.")') thisnam
      write(*,*)
      call terminate
    end if

    ! When coming from sreen->df->dfq->getpemat the bands are set the following way:
    !   occupied: nst1 = sto1-sta1+1, istl1 = sta1, istlu1 = sto1
    !   unoccupied: nst2 = sto2-sta2+1, istl2 = istunocc0+sta2-1, istlu2 = istunocc+sto2-1
    !   and analogously 3 -> unoccupied 4 -> occupied, i.e 34 == 21
    
    ! Gamma point -> Momentum matrix
    gamma: if(tq0) then

      ! Read momentum matrix elements for 12 combination
      call getpmat(ik, vkl0, istl1, istu1, istl2, istu2, .true., trim(pfilnam), p12)

      if(present(p34)) then
        ! Read momentum matrix elements for 34 combination
        call getpmat(ik, vkl0, istl3, istu3, istl4, istu4, .true., trim(pfilnam), p34)
      end if

      ! Consider symmetric gauge w.r.t. coulomb potential
      ! (multiply with v^(1/2))
      ! and normalize w.r.t. KS eigenvalues (no scissors correction!)
      directions: do j = 1, 3

        ! p12
        do i1 = 1, nst1
          do i2 = 1, nst2

            ! epsdfde defaults to 1.0d-8
            if(abs(deou(i1, i2)) .ge. input%xs%epsdfde) then

              ! The square root of 4*pi stems from the square root of the coulomb 
              ! potential
              p12(j, i1, i2) = -p12(j, i1, i2) / deou(i1, i2) * fourpisqt

            else

              p12(j, i1, i2) = zzero

              ! Debug output: Warn if there is a non vanishing occupation difference 
              ! but a vanishing enegy difference in the denominator.
              if((abs(docc12(i1, i2)) .gt. input%groundstate%epsocc)&
                & .and.(input%xs%dbglev .gt. 0)) then

                write(*, '("Warning(", a, "): Divergent energy denominator:&
                  & q-point, k-point, band indices 1-2, delta e12, delta occ:",&
                  & 4i6, 2g18.10)') thisnam, iq, ik, i1 + istl1 - 1,&
                  & i2 + istl2 - 1, deou(i1, i2), docc12(i1,i2)

              end if

            end if

          end do
        end do

        ! p34
        if(present(p34)) then

          do i3 = 1, nst3
            do i4 = 1, nst4

              if(abs(deuo(i3, i4)) .ge. input%xs%epsdfde) then

                p34(j, i3, i4) = - p34(j, i3, i4) / deuo(i3, i4) * fourpisqt

              else

                p34(j, i3, i4) = zzero

                if((abs(docc21(i3, i4)) .gt. input%groundstate%epsocc)&
                  & .and.(input%xs%dbglev .gt. 0)) then

                  write(*, '("Warning(", a, "): Divergent energy denominator:&
                    & q-point, k-point, band indices 3-4:", 4i6, g18.10)')&
                    & thisnam, iq, ik, i4 + istl4 - 1, i3 + istl3 - 1, deuo(i3, i4)

                end if

              end if

            end do
          end do

        end if

      end do directions

    end if gamma

    ! If not G=q=0
    Gqnot0: if(.not. tq0 .or. n .gt. 1) then

      ! For BSE(-kernel) matrix elements are calculated on the fly
      if(tscreen) then

        ! Read plane wave matrix elements form modxs
        m12(:, :, :) = xiou(:, :, :)
        if(present(m34)) m34(:, :, :) = xiuo(:, :, :)

      else

        ! Read matrix elements of plane wave form file
        if(present(m34)) then
          call getemat(iq, ik, .true., trim(efilnam), ngq(iq),&
            & istl1, istu1, istl2, istu2, m12, istl3, istu3,&
            & istl4, istu4, m34)
        else
          call getemat(iq, ik, .true., trim(efilnam), ngq(iq),&
            & istl1, istu1, istl2, istu2, m12)
        end if

      end if

      ! Consider symmetric gauge wrt. Coulomb potential (multiply with v^(1/2))
      ! G not zero
      if(.not. tq0) then

        ! M_ou(G=0,q/=0) = M_ou(G=0,q/=0) * v^1/2(G=0,q/=0)
        m12(:, :, 1) = m12(:, :, 1) *sptclg(1,iq)
        ! M_uo(G=0,q/=0) = M_uo(G=0,q/=0) * v^1/2(G=0,q/=0)
        if(present(m34)) m34(:, :, 1) = m34(:, :, 1) *sptclg(1,iq)

      end if

      ! G+q /= 0
      if(n .gt. 1) then

        forall(igq=2:n)
          m12(:, :, igq) = m12(:, :, igq) *sptclg(igq,iq)
        end forall

        if(present(m34)) then

          forall(igq=2:n)
            m34(:, :, igq) = m34(:, :, igq) *sptclg(igq,iq)
          end forall

        end if

      end if

    end if Gqnot0

  end subroutine getpemat

end module m_getpemat
!EOC
