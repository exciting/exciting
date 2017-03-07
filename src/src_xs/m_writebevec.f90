module m_writebevec

  use modinput
  use modmpi
  use mod_constants, only: h2ev
  use m_getunit
  use m_genfilname
  use m_putgetexcitons

  implicit none

  contains 

    subroutine b_writebevec

      integer(4) :: iqmt, alpha, lambda, nvmax, ncmax, icmin, icmax, ivmin, ivmax 
      real(8) :: rbevec, abevec, en1, en2
      real(8), parameter :: epslat = 1.0d-6
      real(8), allocatable :: rbevec_ksum(:,:), abevec_ksum(:,:)
      integer(4) :: i1, i2, iv, ic, iknr, ikpnr
      integer :: un

      character(256) :: fname, lambdastring

      if(rank == 0) then 

        ! Set defaults if writeexcitons is not specified
        if( .not. associated(input%xs%writeexcitons)) then
          input%xs%writeexcitons => getstructwriteexcitons(emptynode)
        end if

        !====================================================!
        ! Read in data to putgetexcitons module              !
        !====================================================!
        ! Requested iqmt
        iqmt = input%xs%bse%iqmt
        if(iqmt == -1) then
          write(*,*) "Warning(b_writebvec): iqmt=-1, setting it to 1"
          iqmt = 1
        end if
        ! Requested excition index range
        if(input%xs%writeexcitons%selectenergy) then 
          en1=input%xs%writeexcitons%minenergyexcitons
          en2=input%xs%writeexcitons%maxenergyexcitons
          if(input%xs%storeexcitons%useev) then 
            en1=en1/h2ev
            en2=en2/h2ev
          end if
          write(*,*) "m_writebevec: en1, en2", en1, en2
          call get_excitons(iqmt=iqmt, e1=en1, e2=en2)
        else
          i1 = input%xs%writeexcitons%minnumberexcitons
          i2 = input%xs%writeexcitons%maxnumberexcitons
          call get_excitons(iqmt=iqmt, a1=i1, a2=i2)
        end if
        !====================================================!

        !====================================================!
        ! Write out squared modulus of the resonant          !
        ! (anti-resonant) exciton coefficients to text file. !
        !====================================================!

        nvmax = maxval(koulims_(4,:)-koulims_(3,:)+1)
        ncmax = maxval(koulims_(2,:)-koulims_(1,:)+1)
        icmin = minval(koulims_(1,:))
        icmax = maxval(koulims_(2,:))
        ivmin = minval(koulims_(3,:))
        ivmax = maxval(koulims_(4,:))

        ! Loop over eigenvectors
        lambdaloop: do lambda = iex1_, iex2_

          call getunit(un)

          write(lambdastring, '("_LAMBDA",i4.4)') lambda
          fname=trim('BEVEC'//trim(lambdastring)//'.OUT')

          open(unit=un, file=trim(fname), form='formatted', action='write')

          ! Loop over transitions
          do alpha = 1, hamsize_

            ! Get absolute indices form combinded index
            ic = smap_(1,alpha)
            iv = smap_(2,alpha)
            iknr = smap_(3,alpha)

            ikpnr = ikmapikq_(iknr)

            ! Resonant
            rbevec = rvec_(alpha, lambda) * conjg(rvec_(alpha, lambda))
            if(fcoup_) then 
              ! Anti-resonant
              abevec = avec_(alpha, lambda) * conjg(avec_(alpha, lambda))
              write(un, '(2(i8, 3E14.7), 2i8, 2E14.7)')&
                & iknr, vkl0_(:, iknr), ikpnr, vkl_(:,ikpnr),&
                & iv, ic, rbevec, abevec
            else
              write(un, '(2(i8, 3E14.7), 2i8, E14.7)')&
                & iknr, vkl0_(:, iknr), ikpnr, vkl_(:,ikpnr),&
                & iv, ic, rbevec
            end if

          end do

          write(un,*)
          if(fcoup_) then 
            write(un, '("k-point, k-point coordinates, kp-point, kp-point coordinates,&
              & valence band, conduction band, abs(X_alpha)^2, abs(Y_alpha)^2 ")')
          else
            write(un, '("k-point, k-point coordinates, kp-point, kp-point coordinates,&
              & valence band, conduction band, abs(X_alpha)^2 ")')
          end if
          write(un,*)
          write(un, '(i8, " : nr. k-points (total)")') nk_max_
          write(un, '(i8, " : nr. k-points (participating)")') nk_bse_
          write(un, '(i8, " : first (partial) unoccupied band")') iuref_
          write(un, '(i8, " : nr. valence states (maximal) ")') nvmax 
          write(un, '(i8, " : nr. conduction states (maximal)")') ncmax
          write(un,*)

          close (un)

        end do lambdaloop
        !====================================================!

        !====================================================!
        ! Write out squared modulus of the resonant          !
        ! (anti-resonant) exciton coefficients summed over   !
        ! all k-points to text file.                         ! 
        !====================================================!

        allocate(rbevec_ksum(icmin:icmax, ivmin:ivmax))
        if(fcoup_) then 
          allocate(abevec_ksum(icmin:icmax, ivmin:ivmax))
        end if

        do lambda = iex1_, iex2_

          rbevec_ksum=0.0d0
          if(fcoup_) abevec_ksum=0.0d0

          do alpha=1, hamsize_
            ic = smap_(1,alpha)
            iv = smap_(2,alpha)
            rbevec_ksum(ic, iv) = rbevec_ksum(ic, iv) +&
              & rvec_(alpha, lambda)*conjg(rvec_(alpha, lambda))
            if(fcoup_) then 
              abevec_ksum(ic, iv) = abevec_ksum(ic, iv) +&
                & avec_(alpha, lambda)*conjg(avec_(alpha, lambda))
            end if
          end do

          call getunit(un)

          write(lambdastring, '("_LAMBDA",i4.4)') lambda
          fname=trim('BEVEC_KSUM'//trim(lambdastring)//'.OUT')

          open(unit=un, file=trim(fname), form='formatted', action='write')

          do iv=ivmin, ivmax
            do ic=icmin, icmax
              if(fcoup_) then 
                write(un, '(2I8, 2E14.7)') iv, ic, rbevec_ksum(ic, iv), abevec_ksum(ic, iv)
              else
                write(un, '(2I8, E14.7)') iv, ic, rbevec_ksum(ic, iv)
              end if
            end do
          end do

          write(un,*)
          if(fcoup_) then 
            write(un, '("valence band, conduction band, sum_k(abs(X)^2), sum_k(abs(Y)^2)")')
          else
            write(un, '("valence band, conduction band, sum_k(abs(X)^2), sum_k(abs(Y)^2)")')
          end if
          write(un,*)
          write(un, '(i8, " : first (partial) unoccupied band")') iuref_ 
          write(un, '(i8, " : nr. valence states (maximal)")') nvmax
          write(un, '(i8, " : nr. conduction states (maximal)")') ncmax
          write(un,*)

          close(un)

        end do

        deallocate(rbevec_ksum)
        if(fcoup_) then 
          deallocate(abevec_ksum)
        end if
        !====================================================!

        !====================================================!
        ! Legacy output for qmt=0, TDA and fixed bands       !
        !====================================================!
        !---------------------------------------------------------------------------
        ! din: new output file for the bandstructure to be able to post-process it
        !---------------------------------------------------------------------------
        if(iq_ ==1 .and. fcoup_ == .false. .and. fesel_ == .false.) then
          ! Loop over eigenvectors
          do lambda = iex1_, iex2_

            call getunit(un)

            write(fname, '("exciton_evec_",i4.4,".dat")') lambda

            open(unit=un, file=trim(fname), form='formatted', action='write')

            ! nkpt total, nv, iv0, nc, ic0
            write(un,*) "# ", nk_bse_, &
            &                 nvmax, ioref_,  &
            &                 ncmax, iuref_

            ! Loop over transitions
            do alpha = 1, hamsize_

              ! Get absolute indices form combinded index
              ic = smap_(1,alpha)
              iv = smap_(2,alpha)
              iknr = smap_(3,alpha)

              ! Resonant
              rbevec = rvec_(alpha, lambda) * conjg(rvec_(alpha, lambda))
              write(un, '(i8, 4x, 3f10.6, 2i8, g18.10)')&
                & iknr, vkl0_(:, iknr), iv, ic, rbevec

            end do

            close (un)

          end do
        end if
        !====================================================!

        ! Clear read in arrays
        call clear_excitons()

        call barrier

      else

        write(*,*) "m_writebevec: rank", rank, " waiting..."
        call barrier

      end if

    end subroutine b_writebevec

end module m_writebevec
