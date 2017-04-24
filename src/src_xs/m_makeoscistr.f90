module m_makeoscistr

  implicit none

  contains

    subroutine makeoscistr(iqmt, nexc, bevecr, oscstrr, oscstra, bevalre, cmat)
      use modinput
      use modxs, only: unitout, vqlmt, ivgigq, ivgmt
      use modbse, only: hamsize, ensortidx
      use mod_constants
      use m_setup_rmat
      use m_setup_pwmat
      use m_invertzmat
      
      implicit none

      ! I/O
      integer(4), intent(in) :: iqmt, nexc
      complex(8), intent(inout) :: bevecr(:,:)
      complex(8), intent(out) :: oscstrr(nexc,3)
      complex(8), intent(out), optional :: oscstra(nexc, 3)
      real(8), intent(in), optional :: bevalre(:)
      complex(8), intent(inout), optional :: cmat(:, :)

      ! Local
      logical :: useip, useti, usecoup
      real(8) :: t1, t0, ts0, ts1
      complex(8), allocatable :: projmat(:,:), xpy(:,:), rbarmat(:,:)
      integer(4) :: a1, lambda, igqmt, nopt
      integer(4) :: i

      character(*), parameter :: thisname = "makeoscistr"

      ! Use independent particle approximation
      if(input%xs%bse%bsetype == "IP") then
        useip = .true.
      else
        useip = .false.
      end if
      ! Include coupling terms
      usecoup = input%xs%bse%coupling
      ! Use time reversed anti-resonant basis
      useti = input%xs%bse%ti
      
      write(unitout, '("Info(",a,"): Making oscillator strengths.")') trim(thisname)
      write(unitout, '("Info(",a,"): iqmt=", i4)') trim(thisname),  iqmt
      write(unitout, '("Info(",a,"): vqlmt=", 3E11.3)') trim(thisname),  vqlmt(:,iqmt)
      if(useip) then
        write(unitout, '("Info(",a,"): Using IP approximation")') trim(thisname)
      end if

      call timesec(ts0)

      if(iqmt == 1) then 

        write(unitout, '("Info(",a,"): Using zero momentum transfer formalismus.")')&
          & trim(thisname)
        write(unitout, '("Info(",a,"):&
          & Making oscillator strengths using position operator matrix elements.")')&
          & trim(thisname)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! Building position operator matrix elements using momentum matrix elements  !
        ! and transition energies. If on top of GW, renormalize the p mat elements.  !
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        nopt=3
        allocate(projmat(hamsize, nopt))
        projmat=zzero
        call setup_rmat(projmat)

      ! qmt /= 0
      else

        write(unitout, '("Info(",a,"): Making oscillator strengths using&
          & plane wave matrix elements.")') trim(thisname)

        nopt=1
        allocate(projmat(hamsize,nopt))
        ! Generate plane wave matrix elements for G=Gqmt and q=qmt
        igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
        call setup_pwmat(projmat, iqmt, igqmt)

      end if

      ! In case of IP the eigenvectors are not computed, manually
      ! matching indexing to eigenvalue sorting.
      if(useip) then 
        projmat(:,:) = projmat(ensortidx,:)
      end if

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! TDA case: Build resonant oscillator strengths                              !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      if(.not. usecoup) then 
        write(unitout, '("  Building resonant oscillator strengths.")')
        call timesec(t0)
        oscstrr = zzero
        !! Resonant oscillator strengths
        if(.not. useip) then 

          ! qmt=0 case:
          ! t^R_{\lambda,i} = <X_\lambda|\tilde{R}^{i*}> =
          !   \Sum_{a} X^H_{\lambda, a} \tilde{R}^*_{a,i}
          ! qmt/=0 case:
          ! t^R_{\lambda}(G,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
          !   \Sum_{a} X^H_{\lambda, a} \tilde{M}^*_{a}(G,qmt)
          call zgemm('c','n', nexc, nopt, hamsize,&
            & zone, bevecr(1:hamsize,1:nexc), hamsize, conjg(projmat), hamsize,&
            & zzero, oscstrr(1:nexc,1:nopt), nexc)
          call timesec(t1)
          write(unitout, '("    Time needed",f12.3,"s")') t1-t0

        else

          ! IP: X_{\alpha,\lambda} = \delta_{\alpha,\lambda}
          ! qmt=0 case:
          ! t^R_{\lambda,i} = <X_\lambda|\tilde{R}^{i*}> =
          !   \tilde{R}^*_{\lambda,i}
          ! qmt/=0 case:
          ! t^R_{\lambda}(G,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
          !   \tilde{M}^*_{\lambda}(G,qmt)
          oscstrr(1:nexc,1:nopt) = conjg(projmat(1:nexc,1:nopt))

        end if

      end if
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Full case: Build left and right oscillator strengths                       !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      if(usecoup .and. .not. useti) then 

        allocate(rbarmat(2*hamsize,3))

        if(.not. useip) then 
          rbarmat(1:hamsize,1:3) = conjg(projmat)
          rbarmat(hamsize+1:2*hamsize,1:3) = -projmat
        else
          ! In case of IP the eigenvectors are not computed, so 
          ! the order needs to be adapted to match eigenvalue sorting.
          rbarmat(hamsize+1:2*hamsize,1:3) = conjg(projmat)
          do i = 1, hamsize
            rbarmat(i,1:3) = -projmat(hamsize-i+1,:)
          end do
        end if

        deallocate(projmat)

        write(unitout, '("  Building t_r oscillator strengths.")')
        call timesec(t0)

        if(.not. useip) then 
          !! Oscillator strengths from right eigenvectors
          ! t_r_{\lambda,i} = <P_\lambda|\bar{R}_i> =
          !   =  \Sum_{a} P^H_{\lambda,a} \bar{R}_{a,i}
          call zgemm('c','n', nexc, 3, 2*hamsize,&
            & zone, bevecr(1:2*hamsize,1:nexc), 2*hamsize, rbarmat, 2*hamsize,&
            & zzero, oscstrr, nexc)
        else
          ! IP: P_{\alpha,\lambda} = \delta_{\alpha,\lambda}
          ! t_r_{\lambda,i} = <P_\lambda|\bar{R}_i> = \bar{R}_{\lambda,i}
          oscstrr(1:nexc,1:nopt) = rbarmat(1:nexc,1:nopt)
        end if

        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0

        if(.not. useip) then 
          write(unitout, '("  Inverting right EV matrix.")')
          call timesec(t0)
          call zinvert(bevecr)
          call timesec(t1)
          write(unitout, '("    Time needed",f12.3,"s")') t1-t0
        end if

        write(unitout, '("  Building t_l oscillator strengths.")')
        call timesec(t0)

        if(.not. useip) then 
          rbarmat(hamsize+1:2*hamsize,:) = -rbarmat(hamsize+1:2*hamsize,:)
        else
          rbarmat(1:hamsize,:) = -rbarmat(1:hamsize,:)
        end if

        if(.not. useip) then
          !! Oscillator strengths from left eigenvectors
          ! t_l_{\lambda,i} = < P^-1_\lambda | \bar{R}^s_i> =
          !   \Sum_{a} P^{-1}_{\lambda, a} \tilde{R}^s_{a,i}
          call zgemm('n','n', nexc, 3, 2*hamsize,&
            & zone, bevecr(1:nexc,1:2*hamsize), 2*hamsize,&
            & rbarmat, 2*hamsize, zzero, oscstra, nexc)
        else
          ! IP: P^{-1}_{\alpha,\lambda} = \delta_{\alpha,\lambda}
          ! t_l_{\lambda,i} = < P^-1_\lambda | \bar{R}^s_i> =
          !   \tilde{R}^s_{\lambda,i}
          oscstra(1:nexc,1:nopt) = rbarmat(1:nexc, 1:nopt)
        end if

        deallocate(rbarmat)
        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0

      end if
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! TI case: Build oscillator strength                                         !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      if(usecoup .and. useti) then 

        if(.not. useip) then 

          if(.not. present(bevalre) .or. .not. present(cmat)) then
            write(unitout, '("Error(",a,"): Pass eigenvaues and cmat!")') trim(thisname)
            call terminate
          end if

          write(unitout, '("  Building (X+Y) from squared EVP EVs.")')
          call timesec(t0)

          ! Interested in X^+ + Y^+, so we rescale the 
          ! auxiliary eigenvectors Z by the square root of the eigenvalues
          ! so that 
          ! (X+Y)_{a,lambda} = \Sum_{a'} (A-B)^{1/2}_{a,a'} * E^{-1/2}_lambda * Z_{a',lambda}
          do lambda = 1, nexc
            do a1 = 1, hamsize
              bevecr(a1, lambda) = bevecr(a1, lambda) / sqrt(bevalre(lambda))
            end do
          end do
          allocate(xpy(hamsize, nexc))
          call zgemm('N','N',hamsize, nexc, hamsize, zone, cmat, hamsize,&
            & bevecr(:,1:nexc), hamsize, zzero, xpy, hamsize)
          call timesec(t1)
          write(unitout, '("    Time needed",f12.3,"s")') t1-t0

        end if

        write(unitout, '("  Building oscillator strengths for time inverted ar basis.")')
        call timesec(t0)

        if(.not. useip) then
          !! Oscillator strengths
          ! qmt=0
          ! t_{\lambda,i} = < (| X_\lambda>+| Y_\lambda>)| \tilde{R}^{i*}> =
          !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{R}^*_{a,i}
          ! qmt/=0
          ! t_{\lambda}(G,qmt) = < (| X_\lambda>+| Y_\lambda>)| \tilde{M}^*(G,qmt)> =
          !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{M}^*_{a}
          call zgemm('c','n', nexc, nopt, hamsize,&
            & zone, xpy(1:hamsize,1:nexc), hamsize, conjg(projmat), hamsize,&
            & zzero, oscstrr(1:nexc,1:nopt), nexc)
        else
          ! IP= (X+Y)_{\alpha,\lambda} = \delta_{\alpha,\lambda}
          ! qmt=0
          ! t_{\lambda,i} = < (| X_\lambda>+| Y_\lambda>)| \tilde{R}^{i*}> =
          !     \tilde{R}^*_{\lambda,i}
          ! qmt/=0
          ! t_{\lambda}(G,qmt) = < (| X_\lambda>+| Y_\lambda>)| \tilde{M}^*(G,qmt)> =
          !     \tilde{M}^*_{\lambda}
          oscstrr(1:nexc,1:nopt) = conjg(projmat(1:nexc, 1:nopt))
        end if

        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0
        if(allocated(xpy)) deallocate(xpy)
      end if
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      if(allocated(projmat)) deallocate(projmat)
      call timesec(ts1)
      write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0
    end subroutine makeoscistr

    subroutine makeoscistr_dist(iqmt, nexc, dbevecr, doscsr, binfo, bevalre, dcmat)
      use modinput
      use modxs, only: unitout, vqlmt, ivgigq, ivgmt
      use modbse, only: hamsize, ensortidx
      use mod_constants
      use m_setup_rmat
      use m_setup_pwmat
      use m_dzmatmult
      
      implicit none

      ! I/O
      integer(4), intent(in) :: iqmt, nexc
      type(dzmat), intent(inout) :: dbevecr
      type(dzmat), intent(inout) :: doscsr
      type(blacsinfo), intent(in) :: binfo
      real(8), intent(in), optional :: bevalre(:)
      type(dzmat), intent(inout), optional :: dcmat

      ! Local
      logical :: useip, useti, usecoup
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
      ! Use time reversed anti-resonant basis
      useti = input%xs%bse%ti

      write(unitout, '("Info(",a,"):&
        & Making oscillator strengths (distributed).")') trim(thisname)
      write(unitout, '("Info(",a,"): iqmt=", i4)') trim(thisname), iqmt
      write(unitout, '("Info(",a,"): vqlmt=", 3E15.6)') trim(thisname), vqlmt(:,iqmt)
      call timesec(ts0)

      ! qmt = 0
      if(iqmt == 1) then 

        write(unitout, '("Info(",a,"):&
          & Using zero momentum transfer formalismus.")') trim(thisname)
        write(unitout, '("Info(",a,"):&
          & Making oscillator strengths using position operator matrix elements.")')&
          & trim(thisname)

        nopt=3

        if(binfo%isactive) then 

          ! Distributed oscillator strengths
          call new_dzmat(doscsr, nexc, nopt, binfo, rblck=binfo%mblck, cblck=1)

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          ! Building position operator matrix elements using momentum matrix elements  !
          ! and transition energies. If on top of GW, renormalize the p mat elements.  !
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          ! Distributed position operator matrix
          call new_dzmat(dprojmat, hamsize, nopt, binfo, rblck=binfo%mblck, cblck=1)


          ! Build position operator matrix elements
          write(unitout, '("  Building Rmat.")')
          call timesec(t0)

          ! Build R-matrix from P-matrix
          ! \tilde{R}_{a,i} = 
          !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
          !     P^i_{o_a,u_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a})
          call setup_rmat_dist(dprojmat, binfo)

          call timesec(t1)
          write(unitout, '("    Time needed",f12.3,"s")') t1-t0
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
          ! Distributed position operator matrix
          call new_dzmat(dprojmat, hamsize, nopt, binfo, rblck=binfo%mblck, cblck=1)

        else

          write(unitout, '("Info(",a," at rank ", i3,"):&
            & Not active.")')&
            & trim(thisname), mpiglobal%rank

        end if

        ! Build plane wave matrix elements
        write(unitout, '("  Building Pwmat.")')
        call timesec(t0)

        ! Generate plane wave matrix elements for G=Gmt and q=qmt
        igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
        call setup_pwmat_dist(dprojmat, iqmt, igqmt, binfo)

        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      ! qmt 0 or not
      end if

      if(binfo%isactive) then

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! TDA case: Build resonant oscillator strengths                              !
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        if(.not. usecoup) then

          if(binfo%isroot) then 
            write(unitout, '("  Building distributed resonant oscillator strengths.")')
            call timesec(t0)
          end if

          !! Resonant oscillator strengths
          if(.not. useip) then
            ! qmt=0 case:
            ! t^R_{\lambda,i} = <X_\lambda|\tilde{R}^{i*}> =
            !   \Sum_{a} X^H_{\lambda, a} \tilde{R}^*_{a,i}
            ! qmt/=0 case:
            ! t^R_{\lambda}(G,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
            !   \Sum_{a} X^H_{\lambda, a} \tilde{M}^*_{a}(G,qmt)
            dprojmat%za=conjg(dprojmat%za)
            call dzmatmult(dbevecr, dprojmat, doscsr, transa='C', m=nexc)
            ! Deallocating eigenvectors
            call del_dzmat(dbevecr)
          else
            ! IP: X_{\alpha,\lambda} = \delta_{\alpha,\lambda}
            ! qmt=0 case:
            ! t^R_{\lambda,i} = <X_\lambda|\tilde{R}^{i*}> =
            !   \tilde{R}^*_{\lambda,i}
            ! qmt/=0 case:
            ! t^R_{\lambda}(G,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
            !   \tilde{M}^*_{\lambda}(G,qmt)
            doscsr%za(:,:) = conjg(dprojmat%za(:,:))
          end if
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! TI case: Build oscillator strength                                         !
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        if(usecoup .and. useti) then 

          if(.not. useip) then 

            if(.not. present(bevalre) .or. .not. present(dcmat)) then
              write(unitout, '("Error(",a,"): Pass eigenvaues and dcmat!")') trim(thisname)
              call terminate
            end if

            write(unitout, '("  Building (X+Y) from squared EVP EVs.")')
            call timesec(t0)

            ! Interested in X^+ + Y^+, so we rescale the 
            ! auxiliary eigenvectors Z by the square root of the eigenvalues
            ! so that 
            ! (X+Y)_{a,lambda} = \Sum_{a'} (A-B)^{1/2}_{a,a'} * E^{-1/2}_lambda * Z_{a',lambda}
            do j = 1, dbevecr%ncols_loc
              lambda = dbevecr%c2g(j)
              sqrteval = sqrt(bevalre(lambda))
              if(lambda <= nexc) then 
                do i = 1, dbevecr%nrows_loc
                  dbevecr%za(i, j) = dbevecr%za(i, j) / sqrteval
                end do
              end if
            end do
            call new_dzmat(dxpy, hamsize, nexc, binfo)
            call dzmatmult(dcmat, dbevecr, dxpy, n=nexc)

            ! C and Z matrix no longer needed.
            call del_dzmat(dcmat)
            call del_dzmat(dbevecr)

            call timesec(t1)
            write(unitout, '("    Time needed",f12.3,"s")') t1-t0
          end if

          write(unitout, '("  Building distributed oscillator strengths&
            & for time inverted ar basis.")')
          call timesec(t0)

          if(.not. useip) then 
            !! Oscillator strengths
            ! qmt=0
            ! t_{\lambda,i} = < (| X_\lambda>+| Y_\lambda>)| \tilde{R}^{i*}> =
            !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{R}^*_{a,i}
            ! qmt/=0
            ! t_{\lambda}(G,qmt) = < (| X_\lambda>+| Y_\lambda>)| \tilde{M}^*(G,qmt)> =
            !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{M}^*_{a}
            dprojmat%za=conjg(dprojmat%za)
            call dzmatmult(dxpy, dprojmat, doscsr, transa='C', m=nexc) 
            call del_dzmat(dxpy)
          else
            ! IP= (X+Y)_{\alpha,\lambda} = \delta_{\alpha,\lambda}
            ! qmt=0
            ! t_{\lambda,i} = < (| X_\lambda>+| Y_\lambda>)| \tilde{R}^{i*}> =
            !     \tilde{R}^*_{\lambda,i}
            ! qmt/=0
            ! t_{\lambda}(G,qmt) = < (| X_\lambda>+| Y_\lambda>)| \tilde{M}^*(G,qmt)> =
            !     \tilde{M}^*_{\lambda}
            doscsr%za(:,:) = conjg(dprojmat%za(:,:))
          end if

        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        ! Projection matrix no longer needed
        call del_dzmat(dprojmat)

        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0

        call timesec(ts1)
        write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0

      else

        write(unitout, '("Info(",a," at rank ", i3,"):&
          & Not active.")')&
          & trim(thisname), mpiglobal%rank

      ! is active
      end if

    end subroutine makeoscistr_dist

end module m_makeoscistr
