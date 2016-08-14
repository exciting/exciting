! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_getpemat

  implicit none

contains

  subroutine getpemat(iq, ik, pfilnam, efilnam, m12, m34, p12, p34)

    use modinput, only: input
    use mod_constants, only: fourpi, zzero
    use modxs, only: ngq, vkl0, istl1, istl2,&
                   & istl3, istl4, istu1, istu2,&
                   & istu3, istu4, nst1, nst2,&
                   & deou, docc12, deuo, docc21,&
                   & tscreen, xiou, xiuo, sptclg
#ifdef tetra         
    use modtetra
#endif
    use m_getpmat
    use m_getemat

    implicit none

    ! Arguments
    integer, intent(in) :: iq, ik
    character(*), intent(in) :: pfilnam, efilnam
    complex(8), optional, intent(out) :: m12(:, :, :), p34(:, :, :)
    complex(8), optional, intent(out) :: p12(:, :, :), m34(:, :, :)

    ! Local variables
    character(*), parameter :: thisnam = 'getpemat'
    real(8) :: fourpisqt
    integer :: n, igq, j, i1, i2
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

    ! Gamma q-point
    gamma: if(tq0) then

      ! Read momentum matrix elements
      call getpmat(ik, vkl0, istl1, istu1, istl2, istu2, .true., trim(pfilnam), p12)

      if(present(p34))&
        & call getpmat(ik, vkl0, istl3, istu3, istl4, istu4, .true., trim(pfilnam), p34)

      ! Consider symmetric gauge wrt. coulomb potential
      ! (multiply with v^(1/2))
      ! and normalize wrt. ks eigenvalues (no scissors correction!)
      directions: do j = 1, 3

        do i1 = 1, nst1
          do i2 = 1, nst2

            if(abs(deou(i1, i2)) .ge. input%xs%epsdfde) then

              p12(j, i1, i2) = -p12(j, i1, i2) / deou(i1, i2) * fourpisqt

            else

              p12(j, i1, i2) = zzero

              if((abs(docc12(i1, i2)) .gt. input%groundstate%epsocc)&
                & .and.(input%xs%dbglev .gt. 0)) then

                write(*, '("Warning(", a, "): Divergent energy denominator:&
                  & q-point, k-point, band indices 1-2, delta e12, delta occ:",&
                  & 4i6, 2g18.10)') thisnam, iq, ik, i1 + istl1 - 1,&
                  & i2 + istl2 - 1, deou(i1, i2), docc12(i1,i2)
              end if

            end if

            if(present(p34)) then
              
              if(abs(deuo(i2, i1)) .ge. input%xs%epsdfde) then

                p34(j, i2, i1) = - p34 (j, i2, i1) / deuo(i2, i1) * fourpisqt

              else

                p34(j, i2, i1) = zzero

                if((abs(docc21(i2, i1)) .gt. input%groundstate%epsocc)&
                  & .and.(input%xs%dbglev .gt. 0)) then

                  write(*, '("Warning(", a, "): Divergent energy denominator:&
                    & q-point, k-point, band indices 3-4:", 4i6, g18.10)')&
                    & thisnam, iq, ik, i1 + istl1 - 1, i2 + istl2 - 1, deuo(i2, i1)
                end if

               end if

            end if

          end do
        end do

      end do directions

    end if gamma

    ! If not G=q=0
    Gqnot0: if((.not. tq0) .or. (n .gt. 1)) then

      ! For bse(-kernel) matrix elements are calculated on the fly
      if(tscreen) then

        m12(:, :, :) = xiou(:, :, :)
        if(present(m34)) m34(:, :, :) = xiuo(:, :, :)

      else

        ! Read matrix elemets of plane wave
        if(present(m34)) then
          call getemat(iq, ik, .true., trim(efilnam), ngq(iq),&
            & istl1, istu1, istl2, istu2, m12, istl3, istu3,&
            & istl4, istu4, m34)
        else
          call getemat(iq, ik, .true., trim(efilnam), ngq(iq),&
            & istl1, istu1, istl2, istu2, m12)
        end if

      end if

      ! Consider symmetric gauge wrt. coulomb potential (multiply with v^(1/2))
      ! G not zero
      if(.not. tq0) then

        m12(:, :, 1) = m12(:, :, 1) *sptclg(1,iq)
        if(present(m34)) m34(:, :, 1) = m34(:, :, 1) *sptclg(1,iq)

      end if

      ! q not zero
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
