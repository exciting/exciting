module m_makeoscistr

  implicit none

  contains

    subroutine makeoscistr_dist(iqmt, nexc, dexcvec, doscsr, binfo, bevalre, dcmat)
      use modinput
      use modxs, only: unitout, totalqlmt, ivgigq, ivgmt
      use modbse, only: hamsize, ensortidx
      use mod_constants
      use m_setup_dmat
      use m_setup_pwmat
      use m_dzmatmult
      
      implicit none

      ! I/O
      integer(4), intent(in) :: iqmt, nexc
      type(dzmat), intent(inout) :: dexcvec
      type(dzmat), intent(inout) :: doscsr
      type(blacsinfo), intent(in) :: binfo
      real(8), intent(in), optional :: bevalre(:)
      type(dzmat), intent(inout), optional :: dcmat

      ! Local
      logical :: useip, usecoup
      real(8) :: ts0, ts1, t1, t0, sqrteval
      type(dzmat) :: dprojmat, dxpy
      integer(4) :: i,j, lambda, nopt, igqmt

      character(*), parameter :: thisname = "makeoscistr_dist"

      ! Use independent particle approximation
      if(input%xs%bse%bsetype == "IP") then
        useip = .true.
      else
        useip = .false.
      end if
      ! Include coupling terms
      usecoup = input%xs%bse%coupling

      write(unitout, '("Info(",a,"):&
        & Making oscillator strengths (distributed).")') trim(thisname)
      write(unitout, '("Info(",a,"): iqmt=", i4)') trim(thisname), iqmt
      write(unitout, '("Info(",a,"): Qmt=", 3F15.6)') trim(thisname), totalqlmt(:,iqmt)
      if(useip) then
        write(unitout, '("Info(",a,"): Using IP approximation")') trim(thisname)
      end if
      call timesec(ts0)

      ! qmt = 0
      if(iqmt == 1) then 

        write(unitout, '("Info(",a,"):&
          & Using zero momentum transfer formalismus.")') trim(thisname)
        write(unitout, '("Info(",a,"):&
          & Making oscillator strengths using dipole operator matrix elements.")')&
          & trim(thisname)

        nopt=3

        if(binfo%isactive) then 

          ! Distributed oscillator strengths
          call new_dzmat(doscsr, nexc, nopt, binfo, rblck=binfo%mblck, cblck=1)

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          ! Building dipole operator matrix elements using momentum matrix elements    !
          ! and transition energies. If on top of GW, renormalize the p mat elements.  !
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

          ! Distributed dipole operator matrix
          call new_dzmat(dprojmat, hamsize, nopt, binfo, rblck=binfo%mblck, cblck=1)

          ! Build position operator matrix elements
          write(unitout, '("  Building Dipole matrix.")')
          call timesec(t0)

          ! Build D-matrix from P-matrix
          ! \tilde{D}_{a,j} = 
          !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
          !     i*P^j_{o_a,u_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a})
          call setup_dmat_dist(dprojmat, binfo)
          ! Approximated plane wave matrix elements: M -> i*D
          ! Also we use M^* to project the eigenvectors onto in the following
          dprojmat%za=conjg(zi*dprojmat%za)

          call timesec(t1)
          if (input%xs%BSE%outputlevelnumber == 1) &
            & write(unitout, '("    Time needed",f12.3,"s")') t1-t0
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        else

          write(unitout, '("Info(",a," at rank ", i3,"):&
            & Not active.")')&
            & trim(thisname), mpiglobal%rank

        end if

      ! qmt /= 0
      else

        write(unitout, '("Info(",a,"):&
          & Making oscillator strengths using plane wave matrix elements.")')&
          & trim(thisname)

        nopt=1

        if(binfo%isactive) then 

          ! Distributed oscillator strengths
          call new_dzmat(doscsr, nexc, nopt, binfo, rblck=binfo%mblck, cblck=1)

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          ! Building plane wave matrix elements                                        !
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          ! Distributed resonant plane wave matrix
          call new_dzmat(dprojmat, hamsize, nopt, binfo, rblck=binfo%mblck, cblck=1)

        else

          write(unitout, '("Info(",a," at rank ", i3,"):&
            & Not active.")')&
            & trim(thisname), mpiglobal%rank

        end if

        ! Build plane wave matrix elements 
        ! (the pwe generation is parallel over mpiglobal if distribute=true, meaning
        !  also processes that are not on the 2d blacs grid contribute) 
        write(unitout, '("  Building Plane-wave matrix elements.")')
        call timesec(t0)

        ! Generate resonant plane wave matrix elements for G=Gmt and q=qmt
        !   Get corresponding combined G+q index
        igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
        ! Get matrix elements
        ! \tilde{M}_{a} = 
        !   \sqrt{|f_{o_a,k^+_a}-f_{u_a,k^-_a}|} *
        !     M_{u_a o_a k_a-qmt/2}(Gmt,qmt)
        call setup_pwmat_dist(dprojmat, iqmt, igqmt, binfo)
        ! We use M^* to project the eigenvectors onto in the following
        if(binfo%isactive) then 
          dprojmat%za=conjg(dprojmat%za)
        end if

        call timesec(t1)
        if (input%xs%BSE%outputlevelnumber == 1) &
          & write(unitout, '("    Time needed",f12.3,"s")') t1-t0
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      ! qmt 0 or not
      end if

      if(binfo%isactive) then

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! TDA case: Build resonant oscillator strengths                              !
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        if(.not. usecoup) then

          write(unitout, '("  Building distributed resonant oscillator strengths.")')
          call timesec(t0)

          !! Resonant oscillator strengths
          if(.not. useip) then
            ! qmt=0 case:
            ! t^R_{\lambda,j} = <X_\lambda| (i*\tilde{D}^{j})^*> =
            !   -i \Sum_{a} X^H_{\lambda, a} \tilde{D}^*_{a,j}
            ! qmt/=0 case:
            ! t^R_{\lambda}(Gmt,qmt) = <X_\lambda|\tilde{M}^*(Gmt,qmt)> =
            !   \Sum_{a} X^H_{\lambda, a} \tilde{M}^*_{a}(Gmt,qmt)
            call dzmatmult(dexcvec, dprojmat, doscsr, transa='C', m=nexc)
            ! Deallocating eigenvectors
            call del_dzmat(dexcvec)
          else
            ! IP: X_{\alpha,\lambda} = \delta_{\alpha,\lambda}
            ! qmt=0 case:
            ! t^R_{\lambda,j} = <X_\lambda|(i*\tilde{D}^j)^*> =
            !   -i*\tilde{D}^*_{\lambda,j}
            ! qmt/=0 case:
            ! t^R_{\lambda}(Gmt,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
            !   \tilde{M}^*_{\lambda}(G,qmt)
            doscsr%za(:,:) = dprojmat%za(:,:)
          end if
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! TI case: Build oscillator strength                                         !
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        if(usecoup) then 

          if(.not. useip) then 

            ! Input check
            if(.not. present(bevalre) .or. .not. present(dcmat)) then
              write(unitout, '("Error(",a,"):&
                & Pass eigenvaues and dcmat!")') trim(thisname)
              call terminate
            end if

            write(unitout, '("  Building (X+Y) from squared EVP EVs.")')
            call timesec(t0)

            ! Interested in X + Y, so we rescale the 
            ! auxiliary eigenvectors Z by the square root of the eigenvalues
            ! so that 
            ! (X+Y)_{a,lambda} = \Sum_{a'} (A-B)^{1/2}_{a,a'} * E^{-1/2}_lambda * Z_{a',lambda}
            do j = 1, dexcvec%ncols_loc
              lambda = dexcvec%c2g(j)
              sqrteval = sqrt(bevalre(lambda))
              if(lambda <= nexc) then 
                do i = 1, dexcvec%nrows_loc
                  dexcvec%za(i, j) = dexcvec%za(i, j) / sqrteval
                end do
              end if
            end do
            call new_dzmat(dxpy, hamsize, nexc, binfo)
            call dzmatmult(dcmat, dexcvec, dxpy, n=nexc)

            ! C and Z matrix no longer needed.
            call del_dzmat(dcmat)
            call del_dzmat(dexcvec)

            call timesec(t1)
            write(unitout, '("    Time needed",f12.3,"s")') t1-t0
          end if

          write(unitout, '("  Building distributed oscillator strengths&
            & for time reversed ar basis.")')
          call timesec(t0)

          if(.not. useip) then 
            !! Oscillator strengths
            ! qmt=0
            ! t_{\lambda,j} = <X_\lambda + Y_\lambda| (i*\tilde{D}^j)^*> =
            !     -i * \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{D}^*_{a,j}
            ! qmt/=0
            ! t_{\lambda}(Gmt,qmt) = < X_\lambda + Y_\lambda| \tilde{M}^*(Gmt,qmt)> =
            !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{M}^*_{a}(Gmt,qmt)
            call dzmatmult(dxpy, dprojmat, doscsr, transa='C', m=nexc) 
            call del_dzmat(dxpy)
          else
            ! IP= (X+Y)_{\alpha,\lambda} = \delta_{\alpha,\lambda}
            ! qmt=0
            ! t_{\lambda,j} = < X_\lambda + Y_\lambda|(i*\tilde{D}^j)^*> =
            !     -i*\tilde{D}^*_{\lambda,j}
            ! qmt/=0
            ! t_{\lambda}(Gmt,qmt) = < X_\lambda> + Y_\lambda| \tilde{M}^*(Gmt,qmt)> =
            !     \tilde{M}^*_{\lambda}(Gmt,qmt)
            doscsr%za(:,:) = dprojmat%za(:,:)
          end if

        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        ! Projection matrix no longer needed
        call del_dzmat(dprojmat)

        call timesec(t1)
        if (input%xs%BSE%outputlevelnumber == 1) &
          & write(unitout, '("    Time needed",f12.3,"s")') t1-t0

        call timesec(ts1)
        if (input%xs%BSE%outputlevelnumber == 1) &
          & write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0

      else

        write(unitout, '("Info(",a," at rank ", i3,"):&
          & Not active.")')&
          & trim(thisname), mpiglobal%rank

      ! is active
      end if

    end subroutine makeoscistr_dist

end module m_makeoscistr
