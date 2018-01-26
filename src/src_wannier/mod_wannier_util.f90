module mod_wfutil
  use modmain
  use mod_wannier
  use mod_wfint
  use m_plotmat
  implicit none

! module variables

! methods
  contains
    !BOP
    ! !ROUTINE: wfutil_bandstructure
    ! !INTERFACE:
    !
    subroutine wfutil_bandstructure
      ! !USES:
      use m_getunit 
      use FoX_wxml
      ! !INPUT PARAMETERS:
      ! !DESCRIPTION:
      !   Generates Wannier interpolated bandstructure.
      !
      ! !REVISION HISTORY:
      !   Created November 2017 (SeTi)
      !EOP
      !BOC
      integer :: wf_nvp1d, wf_npp1d, iq, iv, lmax, un, is, ia, ias, ist, l
      real(8) :: sum, emin, emax, dk(3), v(3)
      type( k_set) :: tmp_kset
      type( xmlf_t), save :: xf
      character(256) :: fxt, fname, buffer

      real(8), allocatable :: wf_dvp1d(:), wf_vplp1d(:,:), wf_dpp1d(:), bc(:,:,:,:), velo(:,:,:), mass(:,:,:,:)
      real(8), allocatable :: velonum(:), massnum(:), cf(:,:)

      ! generate k-points along path
      wf_nvp1d = size( input%properties%bandstructure%plot1d%path%pointarray)
      wf_npp1d = input%properties%bandstructure%plot1d%path%steps
      allocate( wf_dvp1d( wf_nvp1d))
      allocate( wf_vplp1d( 3, wf_npp1d))
      allocate( wf_dpp1d( wf_npp1d))
      allocate( velonum( wf_npp1d), massnum( wf_npp1d), cf( 3, wf_npp1d))
      call connect( bvec, input%properties%bandstructure%plot1d, wf_nvp1d, wf_npp1d, wf_vplp1d, wf_dvp1d, wf_dpp1d)
      call generate_k_vectors( tmp_kset, bvec, (/1, 1, wf_npp1d/), (/0.d0, 0.d0, 0.d0/), .false.)
      !write(*,*) wf_npp1d, tmp_kset%nkpt
      do iq = 1, tmp_kset%nkpt
        tmp_kset%vkl( :, iq) = wf_vplp1d( :, iq)
        call r3mv( bvec, tmp_kset%vkl( :, iq), tmp_kset%vkc( :, iq))
      end do

      ! interpolate energies
      !call wfint_find_bandgap
      call wfint_init( tmp_kset)
      allocate( velo( 3, wf_nwf, wfint_kset%nkpt))
      allocate( mass( 3, 3, wf_nwf, wfint_kset%nkpt))
      velo = 0.d0
      mass = 0.d0
      call wfint_interpolate_ederiv( velo, mass)
      !write(*,*) wf_npp1d, tmp_kset%nkpt, wfint_kset%nkpt
      !call wfint_interpolate_gwpermat
      !stop

      ! set Fermi energy to zero in output
      fxt = filext
      if( input%properties%wannier%input .eq. "gw") write( filext, '("_GW.OUT")')
      call readfermi
      filext = fxt
      write( fxt, '(".OUT")')
      if( input%properties%wannier%input .eq. "gw") write( fxt, '("_GW.OUT")')
      !wfint_eval = wfint_eval - efermi
      emin = minval( minval( wfint_eval, 1), 1)
      emax = maxval( maxval( wfint_eval, 1), 1)

      ! generate ouput
      if( input%properties%wannier%input .eq. "gw") then
        call xml_OpenFile( "bandstructure-wannier-gw.xml", xf, replace=.true., pretty_print=.true.)
      else
        call xml_OpenFile( "bandstructure-wannier.xml", xf, replace=.true., pretty_print=.true.)
      end if
      call xml_AddXMLPI( xf, "xml-stylesheet", 'href="'//trim(input%xsltpath)//'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
      call xml_NewElement( xf, "bandstructure")
      
      ! interpolate bandcharacter if requested
      if( input%properties%bandstructure%character) then
        lmax = min( 4, input%groundstate%lmaxapw)
        allocate( bc( 0:lmax, natmtot, wf_nwf, wfint_kset%nkpt))
        call wfint_interpolate_bandchar( lmax, bc)
        
        call xml_AddAttribute( xf, "character", "true")
        call xml_NewElement( xf, "title")
        call xml_AddCharacters( xf, trim(input%title))
        call xml_endElement( xf, "title")
        do is = 1, nspecies
          call xml_NewElement( xf, "species")
          call xml_AddAttribute( xf, "name", trim( spname( is)))
          call xml_AddAttribute( xf, "chemicalSymbol", trim( input%structure%speciesarray( is)%species%chemicalSymbol))
          do ia = 1, natoms( is)
            call xml_NewElement( xf, "atom")
            write( buffer, '(5G18.10)') atposc (:, ia, is)
            call xml_AddAttribute( xf, "coord", trim( adjustl( buffer)))
            ias = idxas( ia, is)
            write( fname, '("BAND_WANNIER_S", I2.2, "_A", I4.4)') is, ia
            call getunit( un)
            open( un, file=trim( fname)//trim( fxt), action='write', form='formatted')

            do ist = 1, wf_nwf
              call xml_NewElement( xf, "band")
              do iq = 1, wfint_kset%nkpt
                sum = 0.d0
                do l = 0, lmax
                  sum = sum + bc( l, ias, ist, iq)
                end do
                call xml_NewElement( xf, "point")
                write( buffer, '(5G18.10)') wf_dpp1d( iq)
                call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
                write( buffer, '(5G18.10)') wfint_eval( ist, iq)
                call xml_AddAttribute( xf, "eval", trim( adjustl( buffer)))
                write( buffer, '(5G18.10)') sum
                call xml_AddAttribute( xf, "sum", trim( adjustl( buffer)))
                do l = 0, lmax
                  call xml_NewElement( xf, "bc")
                  write( buffer, *) l
                  call xml_AddAttribute( xf, "l", trim( adjustl( buffer)))
                  write( buffer, '(5G18.10)') bc( l, ias, ist, iq)
                  call xml_AddAttribute( xf, "character", trim( adjustl( buffer)))
                  call xml_endElement( xf, "bc")
                end do
                call xml_endElement( xf, "point")
                write( un, '(2G18.10, 20F12.6)') wf_dpp1d( iq), wfint_eval( ist, iq), sum, bc( :, ias, ist, iq)
              end do
              call xml_endElement( xf, "band")
              write( un, *)
              write( un, *)
            end do
            call xml_endElement( xf, "atom")
            close( un)
          end do
          call xml_endElement( xf, "species")
        end do
        write(*,*)
        write(*, '("Info (wfutil_bandstructure):")')
        write(*,*) "band structure plot written to BAND_WANNIER_Sss_Aaaaa"//trim( fxt)
        write(*, '("	for all species and atoms")')
        deallocate( bc)

      ! without bandcharacter
      else
        call getunit( un)
        open( un, file='BAND_WANNIER'//trim( fxt), action='write', form='formatted')
        call xml_NewElement( xf, "title")
        call xml_AddCharacters( xf, trim( input%title))
        call xml_endElement( xf, "title")
        do ist = 1, wf_nwf
          call xml_NewElement( xf, "band")
          call fderiv( 1, wf_npp1d, wf_dpp1d, wfint_eval( ist, :), velonum, cf)
          call fderiv( 2, wf_npp1d, wf_dpp1d, wfint_eval( ist, :), massnum, cf)
          do iq = 1, wfint_kset%nkpt
            if( iq .eq. wfint_kset%nkpt) then
              dk = wfint_kset%vkc( :, iq) - wfint_kset%vkc( :, iq-1)
            else
              dk = wfint_kset%vkc( :, iq+1) - wfint_kset%vkc( :, iq)
            end if
            dk = dk/norm2( dk)
            call r3mv( mass( :, :, ist, iq), dk, v)
            write( un, '(6G18.10)') wf_dpp1d( iq), wfint_eval( ist, iq), &
                dot_product( velo( :, ist, iq), dk), velonum( iq), &
                dot_product( dk, v), massnum( iq)
            call xml_NewElement( xf, "point")
            write( buffer, '(5G18.10)') wf_dpp1d( iq)
            call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
            write( buffer, '(5G18.10)') wfint_eval( ist, iq)
            call xml_AddAttribute (xf, "eval", trim( adjustl( buffer)))
            call xml_endElement( xf, "point")
          end do
          call xml_endElement( xf, "band")
          write( un, *)
          write( un, *)
        end do
        close( un)
        write(*,*)
        write(*,'("Info (wfutil_bandstructure):")')
        write(*,*) "band structure plot written to BAND_WANNIER"//trim( fxt)
      end if
      write(*,*)
      write(*, '(" Fermi energy is at zero in plot")')

      ! vertex lines
      call getunit( un)
      open( un, file='BANDLINES_WANNIER'//trim( fxt), action='write', form='formatted')
      do iv = 1, wf_nvp1d
        call xml_NewElement( xf, "vertex")
        write( buffer, '(5G18.10)') wf_dvp1d( iv)
        call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
        write( buffer, '(5G18.10)') emax
        call xml_AddAttribute( xf, "upperboundary", trim( adjustl( buffer)))
        write( buffer, '(5G18.10)') emin
        call xml_AddAttribute( xf, "lowerboundary", trim( adjustl( buffer)))
        call xml_AddAttribute( xf, "label", trim( adjustl( input%properties%bandstructure%plot1d%path%pointarray( iv)%point%label))) 
        write( buffer, '(5G18.10)') input%properties%bandstructure%plot1d%path%pointarray( iv)%point%coord
        call xml_AddAttribute( xf, "coord", trim( adjustl( buffer)))
        call xml_endElement( xf, "vertex")
        write( un, '(2G18.10)') wf_dvp1d( iv), emin
        write( un, '(2G18.10)') wf_dvp1d( iv), emax
        write( un, *)
      end do
      close( un)
      write(*,*)
      write(*,*) "vertex location lines written to BANDLINES_WANNIER"//trim( fxt)
      write(*,*)
      call xml_endElement( xf, "bandstructure")
      call xml_close( xf)

      ! bandstructure.dat
      call getunit( un)
      if( input%properties%wannier%input .eq. "gw") then
        open( un, file='bandstructure_wannier_gw.dat', action='write', form='formatted')
      else
        open( un, file='bandstructure_wannier.dat', action='write', form='formatted')
      end if
      write( un, *) "# ", wf_fst, wf_lst, wfint_kset%nkpt
      do ist = 1, wf_nwf
        do iq = 1, wfint_kset%nkpt
          write( un, '(2I6,3F12.6,2G18.10)') ist, iq, wfint_kset%vkl( :, iq), wf_dpp1d( iq), wfint_eval( ist, iq)
        end do
        write( un, *)
      end do
      close( un)

      deallocate( wf_dvp1d, wf_vplp1d, wf_dpp1d)

      return
    end subroutine wfutil_bandstructure

!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfutil_dos
    ! !INTERFACE:
    !
    subroutine wfutil_dos
      ! !USES:
      use m_getunit 
      use FoX_wxml
      ! !INPUT PARAMETERS:
      ! !DESCRIPTION:
      !   Generates Wannier interpolated density of states.
      !
      ! !REVISION HISTORY:
      !   Created November 2017 (SeTi)
      !EOP
      !BOC
      integer :: lmax, lmmax, is, ia, ias, l, m, lm, ie, ntrans, un, n, ist
      integer :: intgrid(3), nsmooth, neffk, nsube
      real(8) :: ewin(2), scissor
      logical :: genpdos, genjdos
      character(256) :: fname, fxt
      character(512) :: buffer
      type (xmlf_t), save :: xf

      real(8), allocatable :: e(:), tdos(:), pdos(:,:,:), jdos(:,:)

      if( minval( input%properties%dos%ngridkint, 1) .le. 0) then
        intgrid = 2*wf_kset%ngridk
        write(*, '(" Warning (wfutil_dos): No or invalid interpolation grid (ngridkint) given: I use twice the original grid.")')
        write(*,'(" given grid: ",3I5)') input%properties%dos%ngridkint
        write(*,'(" used grid:  ",3I5)') intgrid
      else
        intgrid = input%properties%dos%ngridkint
      end if
      
      lmax = min( 3, input%groundstate%lmaxapw)
      lmmax = (lmax+1)**2
      genpdos = input%properties%dos%lmirep
      genjdos = input%properties%dos%jdos
      neffk = input%properties%dos%ngrdos
      nsube = input%properties%dos%nwdos
      nsmooth = input%properties%dos%nsmdos
      scissor = input%properties%dos%scissor
      ewin = input%properties%dos%winddos

      allocate( e( nsube))
      do ie = 1, nsube
        e( ie) = ewin(1) + (ie-0.5d0)*(ewin(2)-ewin(1))/nsube
      end do

      allocate( tdos( nsube))
      if( genpdos) allocate( pdos( nsube, lmmax, natmtot))
      if( genjdos) allocate( jdos( nsube, 0:wf_nwf))

      ! interpolate DOS
      if( genpdos .and. .not. genjdos) then
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             pdos=pdos)
      else if( .not. genpdos .and. genjdos) then
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             jdos=jdos, &
             mtrans=m, &
             ntrans=n)
      else if( genpdos .and. genjdos) then
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             pdos=pdos, &
             jdos=jdos, &
             mtrans=m, &
             ntrans=n)
      else
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor)
      end if

      ! generate output

      ! total DOS
      call getunit( un)
      open( un, file='TDOS_WANNIER'//trim( filext), action='write', form='formatted')
      call xml_OpenFile( "dos_wannier.xml", xf, replace=.True., pretty_print=.True.)
      call xml_NewElement( xf, "dos")
      call xml_NewElement( xf, "title")
      call xml_AddCharacters( xf, trim(input%title))
      call xml_endElement( xf, "title")
      call xml_NewElement( xf, "axis")
      call xml_AddAttribute( xf, "label", 'Energy')
      call xml_AddAttribute( xf, "unit", 'Hartree')
      call xml_endElement( xf, "axis")
      call xml_NewElement( xf, "axis")
      call xml_AddAttribute( xf, "label", 'DOS')
      call xml_AddAttribute( xf, "unit", 'states/Hartree/unit cell')
      call xml_endElement( xf, "axis")

      call xml_NewElement( xf, "totaldos")
      call xml_NewElement( xf, "diagram")
      call xml_AddAttribute( xf, "type", "totaldos")
      write( buffer,*) 1
      call xml_AddAttribute( xf, "nspin", trim( adjustl( buffer)))

      do ie = 1, nsube
        call xml_NewElement (xf, "point")
        write( buffer, '(G18.10)') e( ie)
        call xml_AddAttribute( xf, "e", trim( adjustl( buffer)))
        write( buffer, '(G18.10)') occmax*tdos( ie)
        call xml_AddAttribute( xf, "dos", trim( adjustl( buffer)))
        call xml_endElement( xf, "point")
        write( un, '(2G18.10)') e( ie), occmax*tdos( ie)
      end do
      write( un, *)
      close( un)
      call xml_endElement( xf, "diagram")
      call xml_endElement( xf, "totaldos")

      ! partial and interstitial DOS
      if( genpdos) then
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            call getunit( un)
            write( fname, '("PDOS_WANNIER_S",I2.2,"_A",I4.4)') is, ia
            open( un, file=trim( fname)//trim( filext), action='write', form='formatted')
            call xml_NewElement( xf, "partialdos")
            call xml_AddAttribute( xf, "type", "partial")
            call xml_AddAttribute( xf, "speciessym", trim( adjustl( input%structure%speciesarray(is)%species%chemicalSymbol)))
            write( buffer,*) is
            call xml_AddAttribute( xf, "speciesrn", trim( adjustl( buffer)))
            write( buffer,*) ia
            call xml_AddAttribute( xf, "atom", trim( adjustl( buffer)))
            do l = 0, lmax
              do m = -l, l
                lm = idxlm( l, m)
                call xml_NewElement( xf, "diagram")
                write( buffer,*) 1
                call xml_AddAttribute( xf, "nspin", trim( adjustl( buffer)))
                write( buffer,*) l
                call xml_AddAttribute( xf, "l", trim( adjustl( buffer)))
                write( buffer,*) m
                call xml_AddAttribute( xf, "m", trim( adjustl( buffer)))
                do ie = 1, nsube
                  call xml_NewElement( xf, "point")
                  write( buffer, '(G18.10)') e( ie)
                  call xml_AddAttribute (xf, "e", trim( adjustl( buffer)))
                  write( buffer, '(G18.10)') occmax*pdos( ie, lm, ias)
                  call xml_AddAttribute (xf, "dos", trim( adjustl( buffer)))
                  call xml_endElement( xf, "point")
                  write( un, '(2G18.10)') e( ie), occmax*pdos( ie, lm, ias)
                  tdos( ie) = tdos( ie) - pdos( ie, lm, ias)
                end do
                write( un, *)
                call xml_endElement( xf, "diagram")
              end do
            end do
            close( un)
            call xml_endElement( xf, "partialdos")
          end do
        end do
        call getunit( un)
        open( un, file='IDOS_WANNIER'//trim( filext), action='write', form='formatted')
        call xml_NewElement( xf, "interstitialdos")
        call xml_NewElement( xf, "diagram")
        call xml_AddAttribute( xf, "type", "interstitial")
        write( buffer,*) 1
        call xml_AddAttribute( xf, "nspin", trim( adjustl( buffer)))
        do ie = 1, nsube
          call xml_NewElement( xf, "point")
          write( buffer, '(G18.10)') e( ie)
          call xml_AddAttribute (xf, "e", trim( adjustl( buffer)))
          write( buffer, '(G18.10)') occmax*tdos( ie)
          call xml_AddAttribute (xf, "dos", trim( adjustl( buffer)))
          call xml_endElement( xf, "point")
          write( un, '(2G18.10)') e( ie), occmax*tdos( ie)
        end do
        write( un, *)
        close( un)
        call xml_endElement( xf, "diagram")
        call xml_endElement( xf, "interstitialdos")
        deallocate( pdos)
      end if

      ! joint DOS
      if( genjdos) then
        call getunit( un)
        open( un, file='JDOS_WANNIER'//trim( filext), action='write', form='formatted')
        do ist = 1, n
          do ie = 1, nsube
            if( abs( e( ie)) .gt. 1.d-4) then
              write( un, '(3G18.10)') e( ie), occmax*jdos( ie, ist)/(e( ie)*e( ie))/dble(m*n), occmax*jdos( ie, ist)
            else
              write( un, '(3G18.10)') e( ie), 0.d0, occmax*jdos( ie, ist)
            end if
          end do
          write( un, *)
        end do
        close( un)
        call getunit( un)
        open( un, file='TJDOS_WANNIER'//trim( filext), action='write', form='formatted')
        do ie = 1, nsube
          if( abs( e( ie)) .gt. 1.d-4) then
            write( un, '(3G18.10)') e( ie), occmax*jdos( ie, 0)/( e( ie)*e( ie))/dble(m*n), occmax*jdos( ie, 0)
          else
            write( un, '(3G18.10)') e( ie), 0.d0, occmax*jdos( ie, 0)
          end if
        end do
        write( un, *)
        close( un)
        deallocate( jdos)
      end if        
              
      deallocate( e, tdos)       

      call xml_endElement( xf, "dos")
      call xml_close( xf)

      write(*,*)
      write(*,*) " Info (wfutil_dos):"
      write(*,*)
      write(*,*) "    Total density of states written to TDOS_WANNIER"//trim( filext)
      write(*,*)
      if( genpdos) then
         write(*,*) "    Partial density of states written to PDOS_Sss_Aaaaa"//trim( filext)
         write(*,*) "    for all species and atoms"
         write(*,*)
         write(*,*) "    Interstitial density of states written to IDOS_WANNIER"//trim( filext)
         write(*,*)
      end if
      if( genjdos) then
         write(*,*) "    Joint density of states written to JDOS_WANNIER"//trim( filext)
         write(*,*)
      end if
      write(*,*) "    Fermi energy is at zero in plot"
      write(*,*)
      write(*,*) "    DOS units are states/Hartree/unit cell"
      write(*,*)

      return
    end subroutine wfutil_dos

!--------------------------------------------------------------------------------------

end module mod_wfutil
