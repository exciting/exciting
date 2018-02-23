!BOP
! !ROUTINE: writeexcevec
! !INTERFACE:
subroutine writeexcevec()
! !USES:
  use modxs, only: unitout
  use modinput
  use modmpi
  use mod_constants, only: h2ev
  use m_getunit
  use m_genfilname
  use m_putgetexcitons
! !DESCRIPTION:
!   Reads the binary file \texttt{EXCCOEFF} containing the eigenvector coefficients
!   of the BSE calculation and prints selected coefficients to human readable files.
!   Writes to the folders \texttt{EXCITON\_EVEC}, \texttt{BEVEC} and \texttt{BEVEC\_KSUM}.
!
! !REVISION HISTORY:
!   Created 2016 (Aurich)
!
!EOP
!BOC

  implicit none

  character(*), parameter :: thisname = 'writeexcevec'

  ! Local vars
  integer(4) :: alpha, lambda, nvmax, ncmax, icmin, icmax, ivmin, ivmax 
  real(8) :: rbevec, abevec, en1, en2
  real(8), parameter :: epslat = 1.0d-8
  real(8), allocatable :: rbevec_ksum(:,:), abevec_ksum(:,:)
  integer(4) :: i, i1, i2, iv, ic, ik, ikv, ikc, ivec(3), rcount, acount
  integer(4) :: un
  real(8) :: abscutoffres(2), abscutoffares(2)
  real(8) :: vklv(3), vklc(3)
  real(8), allocatable :: absvec(:)
  integer(4), allocatable :: idxsort(:), idxsort_desc(:)

  integer(4) :: iqmt, iqmti, iqmtf, nqmt, iq1, iq2

  character(256) :: fname, lambdastring
  character(256) :: syscommand, bevecdir, bevecksumdir, excitonevecdir
  character(256) :: tdastring, bsetypestring, scrtypestring

  ! Only rank=0 does the work
  if(mpiglobal%rank == 0) then 

    ! Use Unix system calls to create folders
    bevecdir='BEVEC'
    syscommand = 'test ! -e '//trim(adjustl(bevecdir))//' && mkdir '//trim(adjustl(bevecdir))
    call system(trim(adjustl(syscommand)))
    bevecksumdir='BEVEC_KSUM'
    syscommand = 'test ! -e '//trim(adjustl(bevecksumdir))//' && mkdir '//trim(adjustl(bevecksumdir))
    call system(trim(adjustl(syscommand)))
    excitonevecdir='EXCITON_EVEC'
    syscommand = 'test ! -e '//trim(adjustl(excitonevecdir))//' && mkdir '//trim(adjustl(excitonevecdir))
    call system(trim(adjustl(syscommand)))

    ! Set defaults if writeexcitons is not specified
    if( .not. associated(input%xs%writeexcitons)) then
      input%xs%writeexcitons => getstructwriteexcitons(emptynode)
    end if

    ! Q-point list entries
    !   Use all
    iqmti = 1
    iqmtf = size(input%xs%qpointset%qpoint, 2)
    !   or use range
    if(input%xs%bse%iqmtrange(1) /= -1) then 
      iqmti=input%xs%bse%iqmtrange(1)
      iqmtf=input%xs%bse%iqmtrange(2)
    end if
    nqmt = iqmtf-iqmti+1
    iq1 = 1
    iq2 = nqmt

    ! Info out
    call printline(unitout, "+")
    write(unitout, '("Info(",a,"):", a, i3, a, i3)') trim(thisname),&
      & " Using momentum transfer vectors from list : ", iqmti, " to", iqmtf
    call printline(unitout, "+")

    do iqmt = iqmti+iq1-1, iqmti+iq2-1

      call printline(unitout, "-")
      write(unitout, '("Info(",a,"):", a, i3)') trim(thisname),&
        & " QMT index: ", iqmt

      !====================================================!
      ! Read in data to putgetexcitons module              !
      !====================================================!
      ! Requested excition index range
      !   Use energy selection
      if(input%xs%writeexcitons%selectenergy) then 
        en1=input%xs%writeexcitons%minenergyexcitons
        en2=input%xs%writeexcitons%maxenergyexcitons
        if(input%xs%storeexcitons%useev) then 
          en1=en1/h2ev
          en2=en2/h2ev
        end if
        write(unitout,'("Info(",a,"): Getting excitons in energy range:",2f12.6)')&
          & trim(thisname), en1, en2
        ! Read excition coefficients from binary file
        call get_excitons(iqmt=iqmt, e1=en1, e2=en2)
      !   Use index selection
      else
        i1 = input%xs%writeexcitons%minnumberexcitons
        i2 = input%xs%writeexcitons%maxnumberexcitons
        write(unitout,'("Info(",a,"): Getting excitons in index range:",2i8)')&
          & trim(thisname), i1, i2
        ! Read excition coefficients from binary file
        call get_excitons(iqmt=iqmt, a1=i1, a2=i2)
      end if
      write(unitout,'("Info(",a,"): Excitons read in for range:",2i8)')&
        & trim(thisname), iex1_, iex2_
      !====================================================!

      ! Set maximal number of valence and conduction states over all k-points
      nvmax = maxval(koulims_(4,:)-koulims_(3,:)+1)
      ncmax = maxval(koulims_(2,:)-koulims_(1,:)+1)
      ! Set minimal and maximal valence and conduction state index over all k-points
      icmin = minval(koulims_(1,:))
      icmax = maxval(koulims_(2,:))
      ivmin = minval(koulims_(3,:))
      ivmax = maxval(koulims_(4,:))

      ! Info out
      write(unitout,'("Info(",a,"): nvmax, ncmax",2i8)')&
        & trim(thisname), nvmax, ncmax
      write(unitout,'("Info(",a,"): icmin, icmax",2i8)')&
        & trim(thisname), icmin, icmax
      write(unitout,'("Info(",a,"): ivmin, ivmax",2i8)')&
        & trim(thisname), ivmin, ivmax

      !====================================================!
      ! Selective write of exciton coefficients to file.   !
      !====================================================!

      ! Determine file name from BSE type
      if(fcoup_) then
        tdastring=''
      else
        if(input%xs%bse%chibarq) then
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if
      if(input%xs%bse%bsetype == "IP") then
        tdastring=''
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      ! Printout cutoff for resonant and anti-resonant coefficients
      abscutoffres(1:2) = input%xs%writeexcitons%abscutres
      abscutoffares(1:2) = input%xs%writeexcitons%abscutares

      allocate(absvec(hamsize_))
      allocate(idxsort(hamsize_), idxsort_desc(hamsize_))

      ! For each requested exciton write out the resonant (anti-resonant)
      ! eigen-coefficients sorted by absolute value to text files in folder EXCITON_EVEC
      do lambda = iex1_, iex2_

        write(unitout,'("Info(",a,"): Writing EXEVEC for lambda =",i8)')&
          & trim(thisname), lambda

        call genfilname(dirname=trim(excitonevecdir), basename="EXEVEC",&
          & lambda=lambda, iqmt=iq_,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fname)

        call getunit(un)

        open(unit=un, file=trim(fname), form='formatted', action='write')

        ! Write header
        write(un,'("#",1x,"BSE eigenvector")')
        write(un,'("#")')
        write(un,'("#",1x,"Momentum transfer q_mt:",3f12.7)') vqlmt_
        write(un,'("#",1x,"Size of RR part of the Hamiltonian:",i8)') hamsize_
        write(un,'("#",1x,"Lowest (partially) unoccupied state:",i8)') iuref_
        write(un,'("#")')
        write(un,'("#",1x,"Coefficients are sorted according to thier absolute value")')
        write(un,'("#",2x,"Asolute value cutoff for resonant coefficients:",2E10.3)') abscutoffres
        write(un,'("#",2x,"Asolute value cutoff for anti-resonant coefficients:",2E10.3)') abscutoffares
        write(un,'("#")')
        write(un,'("# Eigenvector number:",i6," with energy/eV: ", f12.7)'), lambda, evals_(lambda)*h2ev
        
        ! Get sorting index for sorting by absolute value
        absvec = abs(rvec_(1:hamsize_,lambda))
        call sortidx(hamsize_, absvec, idxsort)
        do i = 1, hamsize_ 
          idxsort_desc(i) = idxsort(hamsize_-i+1)
        end do

        ! Write resonant coefficients
        write(un,'("# Resonant contribution")')
        write(un,'("#",1x,a,3x,a,7x,a,7x,a,3x,a,33x,a,35x,a,19x,a,15x,a)')&
          & "alpha","lambda", "ic", "iv", "k_c", "k_v",&
          & "|evec|", "real(evec)", "im(evec)"

        ! Note: Resonant basis functions are
        ! \psi^\text{r}_{\alpha,\vec{q}_\text{mt}}(\vec{r},\vec{r}') = 
        !   \phi_{o_\alpha \vec{k}_\alpha+\vec{q}_\text{mt}/2}(\vec{r})
        !   * \phi^*_{u_\alpha \vec{k}_\alpha-\vec{q}_\text{mt}/2}(\vec{r}')

        ! Counter for number of coefficients found in the specified range
        rcount=0

        do i = 1, hamsize_
          ! Get individual indices from combined BSE index
          !   Get combined index from sorted list
          alpha = idxsort_desc(i)
          !   Get corresponding unoccupied state
          ic = smap_(1,alpha)
          !   Get corresponding occupied state
          iv = smap_(2,alpha)
          !   Get reference k-point
          ik = smap_(3,alpha)
          !   Get corresponding k-qmt/2 index and vector
          ikc = ik2ikqmtm_(ik)
          vklc = vkl0_(1:3,ikc)
          !   Get corresponding k+qmt/2 index and vector
          ikv = ik2ikqmtp_(ik)
          vklv = vkl_(1:3,ikv)
          ! Only write out, if the modulus of the coefficient in the required range
          if(absvec(alpha) > abscutoffres(1) .and. absvec(alpha) < abscutoffres(2)) then
            write(un,'(i7,3(2x,i7),2(3f12.7),3(2x,3E25.16))')&
              & alpha, lambda, ic, iv, vklc, vklv,&
              & absvec(alpha), dble(rvec_(alpha,lambda)), aimag(rvec_(alpha, lambda))
            rcount = rcount + 1
          end if
        end do
        write(un,'("# rcount=",i8)'), rcount

        if(fcoup_) then 

          ! Write anti-resonant coefficients
          write(un,'("# Anti-resonant contribution")')
          write(un,'("#",1x,a,3x,a,7x,a,7x,a,3x,a,33x,a,35x,a,19x,a,15x,a)')&
            & "alpha","lambda", "ic", "iv", "k_c", "k_v",&
            & "|evec|", "real(evec)", "im(evec)"

          ! Note: Anti-resonant basis functions are
          ! \psi^\text{a}_{\alpha,\vec{q}_\text{mt}}(\vec{r},\vec{r}') = 
          !   \phi_{u_\alpha -\vec{k}_\alpha+\vec{q}_\text{mt}/2}(\vec{r})
          !   * \phi^*_{o_\alpha -\vec{k}_\alpha-\vec{q}_\text{mt}/2}(\vec{r}')
          ! which is using time reversal symmetry just
          ! \psi^\text{a}_{\alpha,\vec{q}_\text{mt}}(\vec{r},\vec{r}') =
          ! \psi^\text{r}_{\alpha,\vec{q}_\text{mt}}(\vec{r}',\vec{r}) = 
          !   \phi_{o_\alpha \vec{k}_\alpha+\vec{q}_\text{mt}/2}(\vec{r}')
          !   * \phi^*_{u_\alpha \vec{k}_\alpha-\vec{q}_\text{mt}/2}(\vec{r})

          ! Get sorting index for sorting by absolute value
          absvec = abs(avec_(1:hamsize_, lambda))
          call sortidx(hamsize_, absvec, idxsort)
          do i = 1, hamsize_ 
            idxsort_desc(i) = idxsort(hamsize_-i+1)
          end do

          ! Counter for number of coefficients found in the specified range
          acount=0

          do i = 1, hamsize_
            ! Get individual indices from combined BSE index
            !   Get combined index from sorted list
            alpha = idxsort_desc(i)
            !   Get corresponding unoccupied state
            ic = smap_(1,alpha)
            !   Get corresponding occupied state
            iv = smap_(2,alpha)
            !   Get reference k-point
            ik = smap_(3,alpha)
            !   Get corresponding k-qmt/2 index and vector
            ikc = ik2ikqmtm_(ik)
            vklc = vkl0_(1:3,ikc)
            !   Get corresponding k+qmt/2 index and vector
            ikv = ik2ikqmtp_(ik)
            vklv = vkl_(1:3,ikv)

            !vklv = -vklv
            !vklc = -vklc
            !call r3frac(epslat, vklv, ivec)
            !call r3frac(epslat, vklc, ivec)

            ! Only write out, if the modulus of the coefficient in the required range
            if(absvec(alpha) > abscutoffares(1) .and. absvec(alpha) < abscutoffares(2)) then
              write(un,'(i7,3(2x,i7),2(3f12.7),3(2x,3E25.16))')&
                & alpha+hamsize_, lambda, ic, iv, vklc, vklv,&
                & absvec(alpha), dble(avec_(alpha,lambda)), aimag(avec_(alpha, lambda))
              acount=acount+1
            end if
          end do
          write(un,'("# acount=",i8)'), acount

        end if

        close(un)

      end do

      write(unitout,'("Info(",a,"): Writing of EXEVEC finished")')&
        & trim(thisname)

      deallocate(idxsort, idxsort_desc)
      deallocate(absvec)
      !====================================================!

      !====================================================!
      ! Write out squared modulus of the resonant          !
      ! (anti-resonant) exciton coefficients to text file. !
      !====================================================!

      ! Loop over eigenvectors
      lambdaloop: do lambda = iex1_, iex2_

        write(unitout,'("Info(",a,"): Writing BEVEC for lambda =",i8)')&
          & trim(thisname), lambda

        call getunit(un)

        write(lambdastring, '("_LAMBDA",i8.8)') lambda
        fname=trim('BEVEC'//trim(lambdastring)//'.OUT')
        fname=trim(bevecdir)//'/'//trim(fname)

        open(unit=un, file=trim(fname), form='formatted', action='write')

        ! Loop over transitions
        do alpha = 1, hamsize_

          !   Get corresponding unoccupied state
          ic = smap_(1,alpha)
          !   Get corresponding occupied state
          iv = smap_(2,alpha)
          !   Get reference k-point
          ik = smap_(3,alpha)
          !   Get corresponding k-qmt/2 index and vector
          ikc = ik2ikqmtm_(ik)
          vklc = vkl0_(1:3,ikc)
          !   Get corresponding k+qmt/2 index and vector
          ikv = ik2ikqmtp_(ik)
          vklv = vkl_(1:3,ikv)

          ! Resonant
          rbevec = rvec_(alpha, lambda) * conjg(rvec_(alpha, lambda))
          if(fcoup_) then 
            ! Anti-resonant
            abevec = avec_(alpha, lambda) * conjg(avec_(alpha, lambda))
            write(un, '(2(i8, 3E14.7), 2i8, 2E14.7)')&
              & ikv, vklv, ikc, vklc,&
              & iv, ic, rbevec, abevec
          else
            write(un, '(2(i8, 3E14.7), 2i8, 3E14.7)')&
              & ikv, vklv, ikc, vklc,&
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
      write(unitout,'("Info(",a,"): Writing BEVEC finished")')&
        & trim(thisname)
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
        write(unitout,'("Info(",a,"): Writing BEVEC_KSUM for lambda =",i8)')&
          & trim(thisname), lambda

        rbevec_ksum=0.0d0
        if(fcoup_) abevec_ksum=0.0d0

        do alpha=1, hamsize_
          !   Get corresponding unoccupied state
          ic = smap_(1,alpha)
          !   Get corresponding occupied state
          iv = smap_(2,alpha)

          rbevec_ksum(ic, iv) = rbevec_ksum(ic, iv) +&
            & rvec_(alpha, lambda)*conjg(rvec_(alpha, lambda))
          if(fcoup_) then 
            abevec_ksum(ic, iv) = abevec_ksum(ic, iv) +&
              & avec_(alpha, lambda)*conjg(avec_(alpha, lambda))
          end if
        end do

        call getunit(un)

        write(lambdastring, '("_LAMBDA",i8.8)') lambda
        fname=trim('BEVEC_KSUM'//trim(lambdastring)//'.OUT')
        fname=trim(bevecksumdir)//'/'//trim(fname)

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

      ! Clear read in arrays
      call clear_excitons()

    ! iqmt
    end do

    call barrier(callername=trim(thisname))

  else

    write(unitout,'("Info(",a,"): Rank ",i3," is waiting..")')&
      & trim(thisname), mpiglobal%rank

    call barrier(callername=trim(thisname))

  end if

end subroutine writeexcevec
!EOC
