
  subroutine tdcheck
    use modmain, only: reducek, spinsprl, spinpol, version
    use modtddft
    use modtetra
    implicit none
    character(*), parameter :: thisnam = 'tdcheck'
    integer errc, warnc
    character(10) dat, tim

    errc = 0
    warnc = 0

    ! no spin-spirals
    if ( spinsprl ) then
       write(*,*) 'Error('//thisnam//'): not working for spin-spirals'
       errc = errc + 1
    end if

    ! warn for spin polarized calculations
    if ( spinpol ) then
       write(*,*) 'Warning('//thisnam//'): calculation is spin-polarized'
       warnc = warnc + 1
    end if

    ! type of response functions
    if (rsptype.eq.'tord') then
       write(*,'(a,2i8)') 'Error('//thisnam//'): code limitation - only &
            &retarded response functions implemented'
       errc = errc + 1
    end if

    ! tetrahedron method not implemented for analytic continuation
    if (tetra.and.acont) then
       write(*,*) 'Error('//thisnam//'): tetrahedron method does not work &
            &together with analytic continuation'
       errc = errc + 1
    end if

    ! stop on errors
    if ( errc .gt. 0 ) then
       write(*,*) '  Errors occurred - abort'
       call terminate
    end if

    ! warn on warnings
    if ( warnc .gt. 0 ) then
       write(*,*) '  Warnings occurred:', warnc
    end if

  end subroutine tdcheck
