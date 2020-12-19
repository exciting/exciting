module mod_wannier_util
  use modmain
  use mod_wannier
  use mod_wannier_interpolate
  use unit_conversion, only : hartree_to_ev
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
      real(8) :: sum, emin, emax, dk(3), v(3), m(3,3)
      type( k_set) :: tmp_kset
      type( xmlf_t), save :: xf
      character(256) :: fxt, fname, buffer

      real(8), allocatable :: wf_dvp1d(:), wf_vplp1d(:,:), wf_dpp1d(:), bc(:,:,:,:), velo(:,:,:), mass(:,:,:,:), deriv(:,:,:)

      ! generate k-points along path
      wf_nvp1d = size( input%properties%bandstructure%plot1d%path%pointarray)
      wf_npp1d = input%properties%bandstructure%plot1d%path%steps
      allocate( wf_dvp1d( wf_nvp1d))
      allocate( wf_vplp1d( 3, wf_npp1d))
      allocate( wf_dpp1d( wf_npp1d))
      call connect( wf_kset%bvec, wf_nvp1d, wf_npp1d, input%properties%bandstructure%plot1d%path%pointarray, &
        wf_vplp1d, wf_dvp1d, wf_dpp1d)
      call generate_k_vectors( tmp_kset, wf_kset%bvec, (/1, 1, wf_npp1d/), (/0.d0, 0.d0, 0.d0/), .false.)
      !write(*,*) wf_npp1d, tmp_kset%nkpt
      do iq = 1, tmp_kset%nkpt
        tmp_kset%vkl( :, iq) = wf_vplp1d( :, iq)
        call r3mv( wf_kset%bvec, tmp_kset%vkl( :, iq), tmp_kset%vkc( :, iq))
      end do

      ! interpolate energies
      call wfint_init( tmp_kset)
      
      ! interpolate band-derivatives
      if( input%properties%bandstructure%deriv) then
        allocate( deriv( 3, wf_nwf, wfint_kset%nkpt))
        allocate( velo( 3, wf_nwf, wfint_kset%nkpt))
        allocate( mass( 3, 3, wf_nwf, wfint_kset%nkpt))
        call wfint_interpolate_ederiv( velo, mass)
        do ist = 1, wf_nwf
          do iq = 1, wfint_kset%nkpt
            if( iq .eq. 1) then
              dk = wfint_kset%vkc( :, iq+1) - wfint_kset%vkc( :, iq)
            else if( iq .eq. wfint_kset%nkpt) then
              dk = wfint_kset%vkc( :, iq) - wfint_kset%vkc( :, iq-1)
            else
              dk = wfint_kset%vkc( :, iq+1) - wfint_kset%vkc( :, iq-1)
            end if
            dk = dk/norm2( dk)
            call r3mv( mass( :, :, ist, iq), dk, v)
            call r3minv( mass( :, :, ist, iq), m)
            deriv( 1, ist, iq) = dot_product( velo( :, ist, iq), dk)
            deriv( 2, ist, iq) = dot_product( dk, v)
            call r3mv( m, dk, v)
            deriv( 3, ist, iq) = dot_product( dk, v)
          end do
        end do
      end if

      !call wfint_interpolate_eigvec
      !write(*,*) wf_npp1d, tmp_kset%nkpt, wfint_kset%nkpt
      !call wfint_interpolate_gwpermat
      !stop

      ! set Fermi energy to zero in output
      write( fxt, '(".OUT")')
      if( input%properties%wannier%input .eq. "gw") write( fxt, '("_GW.OUT")')
      if( wf_fermizero) wfint_eval = wfint_eval - wf_efermi
      emin = minval( minval( wfint_eval, 1), 1)
      emax = maxval( maxval( wfint_eval, 1), 1)

      ! interpolate bandcharacter if requested
      if( input%properties%bandstructure%character) then
        lmax = min( 4, input%groundstate%lmaxapw)
        allocate( bc( 0:lmax, natmtot, wf_nwf, wfint_kset%nkpt))
        call wfint_interpolate_bandchar( lmax, bc)
      end if

      ! generate ouput
      if( mpiglobal%rank .eq. 0) then
        if( input%properties%wannier%input .eq. "gw") then
          call xml_OpenFile( "bandstructure-wannier-gw.xml", xf, replace=.true., pretty_print=.true.)
        else
          call xml_OpenFile( "bandstructure-wannier.xml", xf, replace=.true., pretty_print=.true.)
        end if
        call xml_AddXMLPI( xf, "xml-stylesheet", 'href="'//trim(input%xsltpath)//'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
        call xml_NewElement( xf, "bandstructure")
        
        ! interpolate bandcharacter if requested
        if( input%properties%bandstructure%character) then
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
              write( buffer, '(5G20.10)') atposc (:, ia, is)
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
                  write( buffer, '(5G20.10)') wf_dpp1d( iq)
                  call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
                  write( buffer, '(5G20.10)') wfint_eval( ist, iq)
                  call xml_AddAttribute( xf, "eval", trim( adjustl( buffer)))
                  write( buffer, '(5G20.10)') sum
                  call xml_AddAttribute( xf, "sum", trim( adjustl( buffer)))
                  if( input%properties%bandstructure%deriv) then
                    write( buffer, '(5G20.10)') deriv( 1, ist, iq) 
                    call xml_AddAttribute( xf, "deriv1", trim( adjustl( buffer)))
                    write( buffer, '(5G20.10)') deriv( 2, ist, iq) 
                    call xml_AddAttribute( xf, "deriv2", trim( adjustl( buffer)))
                    write( buffer, '(5G20.10)') deriv( 3, ist, iq) 
                    call xml_AddAttribute( xf, "mass", trim( adjustl( buffer)))
                  end if
                  do l = 0, lmax
                    call xml_NewElement( xf, "bc")
                    write( buffer, *) l
                    call xml_AddAttribute( xf, "l", trim( adjustl( buffer)))
                    write( buffer, '(5G20.10)') bc( l, ias, ist, iq)
                    call xml_AddAttribute( xf, "character", trim( adjustl( buffer)))
                    call xml_endElement( xf, "bc")
                  end do
                  call xml_endElement( xf, "point")
                  if( input%properties%bandstructure%deriv) then
                    write( un, '(2G20.10, 20F12.6)') wf_dpp1d( iq), wfint_eval( ist, iq), sum, bc( :, ias, ist, iq), deriv( :, ist, iq)
                  else
                    write( un, '(2G20.10, 20F12.6)') wf_dpp1d( iq), wfint_eval( ist, iq), sum, bc( :, ias, ist, iq)
                  end if
                end do
                call xml_endElement( xf, "band")
                write( un, *)
                !write( un, *)
              end do
              call xml_endElement( xf, "atom")
              close( un)
            end do
            call xml_endElement( xf, "species")
          end do
          write(*,*)
          write(*, '("Info (wfutil_bandstructure):")')
          write(*,*) "band structure plot written to BAND_WANNIER_Sss_Aaaaa"//trim( fxt)
          write(*,*) "for all species and atoms"

        ! without bandcharacter
        else
          call getunit( un)
          open( un, file='BAND_WANNIER'//trim( fxt), action='write', form='formatted')
          call xml_NewElement( xf, "title")
          call xml_AddCharacters( xf, trim( input%title))
          call xml_endElement( xf, "title")
          do ist = 1, wf_nwf
            call xml_NewElement( xf, "band")
            do iq = 1, wfint_kset%nkpt
              if( input%properties%bandstructure%deriv) then
                write( un, '(5G20.10)') wf_dpp1d( iq), wfint_eval( ist, iq), deriv( :, ist, iq)
              else
                write( un, '(2G20.10)') wf_dpp1d( iq), wfint_eval( ist, iq)
              end if
              call xml_NewElement( xf, "point")
              write( buffer, '(5G20.10)') wf_dpp1d( iq)
              call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
              write( buffer, '(5G20.10)') wfint_eval( ist, iq)
              call xml_AddAttribute (xf, "eval", trim( adjustl( buffer)))
              if( input%properties%bandstructure%deriv) then
                write( buffer, '(5G20.10)') deriv( 1, ist, iq) 
                call xml_AddAttribute( xf, "deriv1", trim( adjustl( buffer)))
                write( buffer, '(5G20.10)') deriv( 2, ist, iq) 
                call xml_AddAttribute( xf, "deriv2", trim( adjustl( buffer)))
                write( buffer, '(5G20.10)') deriv( 3, ist, iq) 
                call xml_AddAttribute( xf, "mass", trim( adjustl( buffer)))
              end if
              call xml_endElement( xf, "point")
            end do
            call xml_endElement( xf, "band")
            write( un, *)
            !write( un, *)
          end do
          close( un)
          write(*,*)
          write(*,'("Info (wfutil_bandstructure):")')
          write(*,*) "band structure plot written to BAND_WANNIER"//trim( fxt)
        end if
        write(*,'(" Fermi energy is ",g18.10," Hartree.")') wf_efermi
        if( wf_fermizero) then
          write(*,*) "Fermi energy is at zero in plot"
        end if

        ! vertex lines
        call getunit( un)
        open( un, file='BANDLINES_WANNIER'//trim( fxt), action='write', form='formatted')
        do iv = 1, wf_nvp1d
          call xml_NewElement( xf, "vertex")
          write( buffer, '(5G20.10)') wf_dvp1d( iv)
          call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
          write( buffer, '(5G20.10)') emax
          call xml_AddAttribute( xf, "upperboundary", trim( adjustl( buffer)))
          write( buffer, '(5G20.10)') emin
          call xml_AddAttribute( xf, "lowerboundary", trim( adjustl( buffer)))
          call xml_AddAttribute( xf, "label", trim( adjustl( input%properties%bandstructure%plot1d%path%pointarray( iv)%point%label))) 
          write( buffer, '(5G20.10)') input%properties%bandstructure%plot1d%path%pointarray( iv)%point%coord
          call xml_AddAttribute( xf, "coord", trim( adjustl( buffer)))
          call xml_endElement( xf, "vertex")
          write( un, '(2G20.10)') wf_dvp1d( iv), emin
          write( un, '(2G20.10)') wf_dvp1d( iv), emax
          write( un, *)
        end do
        close( un)
        write(*,*) "Vertex location lines written to BANDLINES_WANNIER"//trim( fxt)
        write(*,*)
        call xml_endElement( xf, "bandstructure")
        call xml_close( xf)

        ! bandstructure.dat
        call getunit( un)
        if( input%properties%wannier%input .eq. "gw") then
          open( un, file='bandstructure_wannier_gw.dat', action='write', form='formatted')
        else if( input%properties%wannier%input .eq. "qsgw") then
          open( un, file='bandstructure_wannier_qsgw.dat', action='write', form='formatted')
        else
          open( un, file='bandstructure_wannier.dat', action='write', form='formatted')
        end if
        write( un, *) "# ", wf_fst, wf_fst+wf_nwf-1, wfint_kset%nkpt
        do ist = 1, wf_nwf
          do iq = 1, wfint_kset%nkpt
            write( un, '(2I6,3F12.6,2G20.10)') ist+wf_fst-1, iq, wfint_kset%vkl( :, iq), wf_dpp1d( iq), wfint_eval( ist, iq)
          end do
          write( un, *)
        end do
        close( un)
      end if

      if( allocated( bc)) deallocate( bc)
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
      use mod_optkgrid, only: getoptkgrid
      ! !INPUT PARAMETERS:
      ! !DESCRIPTION:
      !   Generates Wannier interpolated density of states.
      !
      ! !REVISION HISTORY:
      !   Created November 2017 (SeTi)
      !EOP
      !BOC
      integer :: lmax, lmmax, is, ia, ias, l, m, lm, ie, un, ist, ntrans, mtrans, nsm
      integer :: intgrid(3), nsmooth, neffk, nsube
      real(8) :: ewin(2), scissor, rad, opt, ropt, v(3), vbz
      logical :: genpdos, genjdos
      character(64) :: inttype
      character(256) :: fname
      character(512) :: buffer
      type (xmlf_t), save :: xf

      real(8), allocatable :: e(:), tdos(:), pdos(:,:,:), jdos(:,:)

      rad = 1.d-2

      ! determine interpolation grid
      if( input%properties%dos%linkpt .gt. 1.d-16) then
        call getoptkgrid( 1.d-4, wf_kset%bvec, intgrid, opt, ropt)
        call r3cross( wf_kset%bvec(:,2), wf_kset%bvec(:,3), v)
        vbz = dot_product( wf_kset%bvec(:,1), v)
        ropt = product( dble( intgrid))
        intgrid = nint( intgrid*input%properties%dos%linkpt/(ropt/vbz)**(1.d0/3.d0))
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write(*, '(" Info (wfutil_dos): Interpolation grid automatically determined from linkpt = ",f13.6)') input%properties%dos%linkpt
          write(*,'(" used grid:  ",3I5)') intgrid
        end if
      else if( minval( input%properties%dos%ngridkint, 1) .gt. 0) then
        intgrid = input%properties%dos%ngridkint
      else
        call getoptkgrid( 1.d-4, wf_kset%bvec, intgrid, opt, ropt)
        call r3cross( wf_kset%bvec(:,2), wf_kset%bvec(:,3), v)
        vbz = dot_product( wf_kset%bvec(:,1), v)
        ropt = product( dble( intgrid))
        intgrid = nint( intgrid*10.d0/(ropt/vbz)**(1.d0/3.d0))
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write(*, '(" Warning (wfutil_dos): No or invalid interpolation grid (ngridkint) given. Interpolation grid automatically determined.")')
          write(*,'(" given grid: ",3I5)') input%properties%dos%ngridkint
          write(*,'(" used grid:  ",3I5)') intgrid
        end if
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
      inttype = trim( input%properties%dos%inttype)
      nsm = input%properties%dos%nsmdos
      ! smoothing is already done in old integration routine
      if( input%properties%dos%inttype /= 'tetra') nsm = 0

      allocate( e( nsube))
      do ie = 1, nsube
        e( ie) = ewin(1) + dble( ie-1)*(ewin(2)-ewin(1))/(nsube-1)
      end do

      allocate( tdos( nsube))
      if( genpdos) allocate( pdos( nsube, lmmax, natmtot))
      if( genjdos) allocate( jdos( nsube, 0:wf_nwf))

      ! interpolate DOS
      if( genpdos .and. .not. genjdos) then
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             pdos=pdos, &
             inttype=inttype, &
             lonly=input%properties%dos%lonly)
      else if( .not. genpdos .and. genjdos) then
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             jdos=jdos, &
             inttype=inttype, &
             ntrans=ntrans, mtrans=mtrans)
      else if( genpdos .and. genjdos) then
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             pdos=pdos, &
             lonly=input%properties%dos%lonly, &
             jdos=jdos, &
             inttype=inttype, &
             ntrans=ntrans, mtrans=mtrans)
      else
        call wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, &
             scissor=scissor, &
             inttype=inttype)
      end if

      if( mpiglobal%rank .eq. 0) then

        ! generate output
        ! total DOS
        if( nsm .gt. 0) call fsmooth( nsm, nsube, 1, tdos)
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
          write( buffer, '(G20.10E3)') e( ie)
          call xml_AddAttribute( xf, "e", trim( adjustl( buffer)))
          write( buffer, '(G20.10E3)') occmax*tdos( ie)
          call xml_AddAttribute( xf, "dos", trim( adjustl( buffer)))
          call xml_endElement( xf, "point")
          write( un, '(2G20.10E3)') e( ie), occmax*tdos( ie)
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
                if( input%properties%dos%lonly) then
                  if( nsm .gt. 0) call fsmooth( nsm, nsube, 1, pdos( :, l+1, ias))
                  call xml_NewElement( xf, "diagram")
                  write( buffer,*) 1
                  call xml_AddAttribute( xf, "nspin", trim( adjustl( buffer)))
                  write( buffer,*) l
                  call xml_AddAttribute( xf, "l", trim( adjustl( buffer)))
                  do ie = 1, nsube
                    call xml_NewElement( xf, "point")
                    write( buffer, '(G20.10E3)') e( ie)
                    call xml_AddAttribute (xf, "e", trim( adjustl( buffer)))
                    write( buffer, '(G20.10E3)') occmax*pdos( ie, l+1, ias)
                    call xml_AddAttribute (xf, "dos", trim( adjustl( buffer)))
                    call xml_endElement( xf, "point")
                    write( un, '(2G20.10E3)') e( ie), occmax*pdos( ie, l+1, ias)
                    tdos( ie) = tdos( ie) - pdos( ie, l+1, ias)
                  end do
                  write( un, *)
                  call xml_endElement( xf, "diagram")
                else
                  do m = -l, l
                    lm = idxlm( l, m)
                    if( nsm .gt. 0) call fsmooth( nsm, nsube, 1, pdos( :, lm, ias))
                    call xml_NewElement( xf, "diagram")
                    write( buffer,*) 1
                    call xml_AddAttribute( xf, "nspin", trim( adjustl( buffer)))
                    write( buffer,*) l
                    call xml_AddAttribute( xf, "l", trim( adjustl( buffer)))
                    write( buffer,*) m
                    call xml_AddAttribute( xf, "m", trim( adjustl( buffer)))
                    do ie = 1, nsube
                      call xml_NewElement( xf, "point")
                      write( buffer, '(G20.10E3)') e( ie)
                      call xml_AddAttribute (xf, "e", trim( adjustl( buffer)))
                      write( buffer, '(G20.10E3)') occmax*pdos( ie, lm, ias)
                      call xml_AddAttribute (xf, "dos", trim( adjustl( buffer)))
                      call xml_endElement( xf, "point")
                      write( un, '(2G20.10E3)') e( ie), occmax*pdos( ie, lm, ias)
                      tdos( ie) = tdos( ie) - pdos( ie, lm, ias)
                    end do
                    write( un, *)
                    call xml_endElement( xf, "diagram")
                  end do
                end if
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
            write( buffer, '(G20.10E3)') e( ie)
            call xml_AddAttribute (xf, "e", trim( adjustl( buffer)))
            write( buffer, '(G20.10E3)') occmax*tdos( ie)
            call xml_AddAttribute (xf, "dos", trim( adjustl( buffer)))
            call xml_endElement( xf, "point")
            write( un, '(2G20.10E3)') e( ie), occmax*tdos( ie)
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
          do ist = 1, ntrans
            if( nsm .gt. 0) call fsmooth( nsm, nsube, 1, jdos( :, ist))
            do ie = 1, nsube
              if( abs( e( ie)) .gt. 1.d-4) then
                write( un, '(3G20.10E3)') e( ie), jdos( ie, ist)/(e( ie)*e( ie))/dble( ntrans*mtrans), jdos( ie, ist)
              else
                write( un, '(3G20.10E3)') e( ie), 0.d0, jdos( ie, ist)
              end if
            end do
            write( un, *)
          end do
          close( un)
          if( nsm .gt. 0) call fsmooth( nsm, nsube, 1, jdos( :, 0))
          call getunit( un)
          open( un, file='TJDOS_WANNIER'//trim( filext), action='write', form='formatted')
          do ie = 1, nsube
            if( abs( e( ie)) .gt. 1.d-4) then
              write( un, '(3G20.10E3)') e( ie), jdos( ie, 0)/( e( ie)*e( ie))/dble( ntrans*mtrans), jdos( ie, 0)
            else
              write( un, '(3G20.10E3)') e( ie), 0.d0, jdos( ie, 0)
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
        write(*, '("Info (wfutil_dos):")')
        write(*,*) "Total density of states written to TDOS_WANNIER"//trim( filext)
        if( genpdos) then
           write(*,*) "Partial density of states written to PDOS_Sss_Aaaaa"//trim( filext)
           write(*,*) "for all species and atoms"
           write(*,*) "Interstitial density of states written to IDOS_WANNIER"//trim( filext)
        end if
        if( genjdos) then
           write(*,*) "Joint density of states written to JDOS_WANNIER"//trim( filext)
        end if
        write(*,'(" Fermi energy is ",g18.10," Hartree.")') wfint_efermi
        if( wf_fermizero) then
          write(*,*) "Fermi energy is at zero in plot"
        end if
        write(*,*) "DOS units are states/Hartree/unit cell"
        write(*,*)
      end if

      if( allocated( e)) deallocate( e)
      if( allocated( tdos)) deallocate( tdos)
      if( allocated( pdos)) deallocate( pdos)
      if( allocated( jdos)) deallocate( jdos)

      return
    end subroutine wfutil_dos

!--------------------------------------------------------------------------------------
    
    subroutine wfutil_find_bandgap
      use mod_eigenvalue_occupancy, only: occmax
      use mod_charge_and_moment, only: chgval
      use mod_optkgrid, only: getoptkgrid
      use m_getunit

      integer :: nvm, iqvbm, iqcbm, ndiv(3), degvbm, degcbm, ist, n1, n2, un, nkpt
      real(8) :: vvbm(3), vcbm(3)
      real(8) :: v(3), m(3,3), eps_deg, rad, opt, ropt
      type( k_set) :: tmp_kset
      logical :: findvbm, findcbm

      integer, allocatable :: bndlist(:)
      real(8), allocatable :: energy(:), velo(:,:,:), mass(:,:,:,:), kptlist(:,:)
      logical, allocatable :: extlist(:)

      eps_deg = 1.d-4
      rad = 2.d-2

      findvbm = .false.
      findcbm = .false.
      nvm = nint( chgval/occmax)

      if( associated( input%properties%wanniergap)) then
        if( input%properties%wanniergap%auto) then
          if( (wf_fst .ne. 1) .and. (wf_fst .le. nvm)) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Warning (wfutil_find_bandgap): The lowest wannierized band is ",I3,". All bands below are considered to be fully occupied.")') wf_fst
            end if
          end if
          if( wf_fst .gt. nvm) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Warning (wfutil_find_bandgap): No valence bands have been wannierized. Cannot find VBM.")')
            end if
          else
            findvbm = .true.
          end if
          if( (wf_lst .le. nvm)) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Warning (wfutil_find_bandgap): No conduction bands have been wannierized. Cannot find CBM.")')
            end if
          else
            findcbm = .true.
          end if

          nvm = nvm - wf_fst + 1

          if( minval( input%properties%wanniergap%ngridkint) .gt. 0) then
            ndiv = input%properties%wanniergap%ngridkint
          else
            call getoptkgrid( rad, wf_kset%bvec, ndiv, opt, ropt)
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Warning (wfutil_find_bandgap): No interpolation grid given. I use ",3i4,".")') ndiv
            end if
          end if

          call generate_k_vectors( tmp_kset, wf_kset%bvec, ndiv, (/0.d0, 0.d0, 0.d0/), .true., uselibzint=.false.)
          call wfint_init( tmp_kset)

          if( findvbm .and. findcbm) then
            if( maxval( wfint_eval( nvm, :)) .gt. minval( wfint_eval( nvm+1, :))) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wfutil_find_bandgap): I think your system is metalic. No gap can be found.")')
              end if
              return
            end if
          end if

          if( findvbm) then
            iqvbm = maxloc( wfint_eval( nvm, :), 1)
            vvbm = wfint_kset%vkl( :, iqvbm)
          end if

          if( findcbm) then
            iqcbm = minloc( wfint_eval( nvm+1, :), 1)
            vcbm = wfint_kset%vkl( :, iqcbm)
          end if
        end if
      end if

      nkpt = size( input%properties%wanniergap%pointbandarray)
      if( findvbm) nkpt = nkpt + 1
      if( findcbm) nkpt = nkpt + 1
      allocate( kptlist( 3, nkpt))
      allocate( extlist( nkpt), bndlist( nkpt))
      n1 = 0
      extlist = .false.
      if( findvbm) then
        n1 = n1 + 1
        kptlist( :, n1) = vvbm
        bndlist( n1) = nvm
        extlist( n1) = .true.
      end if
      if( findcbm) then
        n1 = n1 + 1
        kptlist( :, n1) = vcbm
        bndlist( n1) = nvm + 1
        extlist( n1) = .true.
      end if
      do ist = 1, size( input%properties%wanniergap%pointbandarray)
        n1 = n1 + 1
        kptlist( :, n1) = input%properties%wanniergap%pointbandarray( ist)%pointband%vkl
        bndlist( n1) = input%properties%wanniergap%pointbandarray( ist)%pointband%band
        extlist( n1) = input%properties%wanniergap%pointbandarray( ist)%pointband%extremal
      end do

      call generate_k_vectors( tmp_kset, wf_kset%bvec, (/1, 1, 1/), (/0.d0, 0.d0, 0.d0/), .false.)
      allocate( energy( nkpt))
      allocate( velo( 3, wf_nwf, nkpt))
      allocate( mass( 3, 3, wf_nwf, nkpt))

      do ist = 1, nkpt
        if( extlist( ist)) call wfutil_energy_extremal( kptlist( :, ist), bndlist( ist), energy( ist), v)
      end do

      call getunit( un)
      open( un, file='GAP_WANNIER.OUT', action='write', form='formatted')

      call generate_k_vectors( tmp_kset, wf_kset%bvec, (/1, 1, nkpt/), (/0.d0, 0.d0, 0.d0/), .false.)
      do ist = 1, nkpt
        tmp_kset%vkl( :, ist) = kptlist( :, ist)
        call r3mv( wf_kset%bvec, tmp_kset%vkl( :, ist), tmp_kset%vkc( :, ist))
      end do
      call wfint_init( tmp_kset)
      call wfint_interpolate_ederiv( velo, mass)

      nvm = 0
      if( findvbm) then
        nvm = nvm + 1
        degvbm = 0
        do ist = 1, wf_nwf
          if( abs( wfint_eval( ist, nvm) - wfint_eval( bndlist( nvm), nvm)) .lt. eps_deg) degvbm = degvbm + 1
        end do
        wfint_evbm = wfint_eval( bndlist( nvm), nvm)

        write( un, *) "VALENCE-BAND MAXIMUM"
        write( un, '(" position (lattice):      ",3f16.6)') wfint_kset%vkl( :, nvm)
        write( un, '("          (cartesian):    ",3f16.6)') wfint_kset%vkc( :, nvm)
        write( un, '(" energy (Hartree):        ",f16.6)') wfint_evbm
        write( un, '("        (eV):             ",f16.6)') wfint_evbm*hartree_to_ev
        write( un, '(" degeneracy:                 ",i6)') degvbm
        do n1 = 1, degvbm
          call r3minv( mass( :, :, bndlist( nvm)-n1+1, nvm), m)
          write( un, '(" band gradient  ",i2,":        ",3g16.6)') n1, velo( :, bndlist( nvm)-n1+1, nvm)
          write( un, '(" effective mass ",i2,":       ",3f16.6)') n1, m( 1, :)
          write( un, '("                          ",3f16.6)') m( 2, :)
          write( un, '("                          ",3f16.6)') m( 3, :)
          write( un, *)
        end do
      end if

      if( findcbm) then
        nvm = nvm + 1
        degcbm = 0
        do ist = 1, wf_nwf
          if( abs( wfint_eval( ist, nvm) - wfint_eval( bndlist( nvm), nvm)) .lt. eps_deg) degcbm = degcbm + 1
        end do
        wfint_ecbm = wfint_eval( bndlist( nvm), nvm)

        write( un, *) "CONDUCTION-BAND MINIMUM"
        write( un, '(" position (lattice):      ",3f16.6)') wfint_kset%vkl( :, nvm)
        write( un, '("          (cartesian):    ",3f16.6)') wfint_kset%vkc( :, nvm)
        write( un, '(" energy (Hartree):        ",f16.6)') wfint_ecbm
        write( un, '("        (eV):             ",f16.6)') wfint_ecbm*hartree_to_ev
        write( un, '(" degeneracy:                 ",i6)') degcbm
        do n1 = 1, degcbm
          call r3minv( mass( :, :, bndlist( nvm)-n1+1, nvm), m)
          write( un, '(" band gradient  ",i2,":        ",3g16.6)') n1, velo( :, bndlist( nvm)-n1+1, nvm)
          write( un, '(" effective mass ",i2,":       ",3f16.6)') n1, m( 1, :)
          write( un, '("                          ",3f16.6)') m( 2, :)
          write( un, '("                          ",3f16.6)') m( 3, :)
          write( un, *)
        end do
      end if

      if( findvbm .and. findcbm) then
        write( un, '(" gap (Hartree):           ",f16.6)') wfint_ecbm - wfint_evbm
        write( un, '("     (eV):                ",f16.6)') (wfint_ecbm - wfint_evbm)*hartree_to_ev
        write( un, *)
      end if

      n2 = nvm + 1
      do nvm = n2, nkpt
        degcbm = 0
        do ist = 1, wf_nwf
          if( abs( wfint_eval( ist, nvm) - wfint_eval( bndlist( nvm), nvm)) .lt. eps_deg) degcbm = degcbm + 1
        end do

        write( un, '(" POINT ",i2)') nvm - n2 + 1
        write( un, '(" band index:                 ",i6)') bndlist( nvm)
        write( un, '(" position (lattice):      ",3f16.6)') wfint_kset%vkl( :, nvm)
        write( un, '("          (cartesian):    ",3f16.6)') wfint_kset%vkc( :, nvm)
        write( un, '(" energy (Hartree):        ",f16.6)') wfint_eval( bndlist( nvm), nvm)
        write( un, '("        (eV):             ",f16.6)') wfint_eval( bndlist( nvm), nvm)*hartree_to_ev
        write( un, '(" degeneracy:                 ",i6)') degcbm
        do n1 = 1, degcbm
          call r3minv( mass( :, :, bndlist( nvm)-n1+1, nvm), m)
          write( un, '(" band gradient  ",i2,":        ",3g16.6)') n1, velo( :, bndlist( nvm)-n1+1, nvm)
          write( un, '(" effective mass ",i2,":       ",3f16.6)') n1, m( 1, :)
          write( un, '("                          ",3f16.6)') m( 2, :)
          write( un, '("                          ",3f16.6)') m( 3, :)
          write( un, *)
        end do
      end do

      close( un)

      write(*,'(" Infos on bandgap written to GAP_WANNIER.OUT.")')

      return
    end subroutine wfutil_find_bandgap

    subroutine wfutil_energy_extremal( vkl, ist, eval, velo)
      real(8), intent( inout)  :: vkl(3)
      integer, intent( in)     :: ist
      real(8), intent( out)    :: eval, velo(3)

      integer, parameter :: maxit = 100
      real(8), parameter :: eps   = 1.d-6

      integer :: i, iv(3)
      real(8) :: s, bi(3,3), vt(3), v( 3, wf_nwf, 1), m( 3, 3, wf_nwf, 1)
      type( k_set) :: kset

      call generate_k_vectors( kset, wf_kset%bvec, (/1, 1, 1/), (/0.d0, 0.d0, 0.d0/), .false.)
      call r3minv( kset%bvec, bi)

      velo = (/1.d0, 1.d0, 1.d0/)
      i = 0
      do while( (norm2( velo) .gt. eps) .and. (i .lt. maxit))
        call r3ws( 1.d-2*eps, kset%bvec, vkl, iv)
        i = i + 1
        kset%vkl( :, 1) = vkl
        call r3mv( kset%bvec, kset%vkl( :, 1), kset%vkc( :, 1))
        call wfint_init( kset)
        call wfint_interpolate_ederiv( v, m)
        velo = v( :, ist, 1)
        call r3mv( m( :, :, ist, 1), velo, vt)
        s = dot_product( velo, velo)/dot_product( velo, vt)
        call r3mv( bi, velo, vt)
        vkl = vkl - s*vt
        eval = wfint_eval( ist, 1)
        !write(*,'(i,3f13.6,f23.16,g20.6)') i, vkl, eval, norm2( velo)
      end do

      call delete_k_vectors( kset)
      return
    end subroutine wfutil_energy_extremal
!--------------------------------------------------------------------------------------

    subroutine wfutil_plot( fst, lst, cell)
      use modplotlabels
      use mod_rgrid
      use mod_xsf_format
      use mod_cube_format
      use m_wsweight
      use mod_lattice, only: omega, ainv
      implicit none
      ! input/output
      integer, intent(in) :: fst, lst, cell(3)
      ! local variables
      integer :: ip, np, nv, ik, iknr, ist, is, ia, ias, maxnpt
      character(80) :: fname
      real(8) :: cellc(3), s, phi, v0(3), v1(3), v2(3), v3(3), rrange( 2, fst:lst), irange( 2, fst:lst) 
      complex(8) :: phase
      ! allocatable arrays
      real(8), allocatable :: dist(:)
  !    complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), dmatk(:,:,:,:,:), dmatr(:,:,:,:)
      complex(8), allocatable :: evecfv(:,:), apwalm(:,:,:,:), evec(:)
      complex(8), allocatable :: wfmt(:,:,:,:), wfir(:,:)
      complex(8), allocatable :: zdata(:), zdatatot(:,:)
      !
      type(rgrid), allocatable :: grid(:)
  
      ! initialise universal variables
      input%groundstate%lradstep = 1
      ! WARNING: make sure that the correct state is read and correct global k-point arrays are definded
      !call init0
      !call init1    
      call wfhelp_genradfun
  
      if ((fst<1) .or. (lst>wf_nwf)) then
        write(*,*)
        write(*, '("Error (wfutil_plot): state out of range : ", I8)') ist
        stop
      end if
  
      call r3mv( input%structure%crystal%basevect, dble( cell), cellc)
      allocate( grid( fst:lst))
  
      ! generate plotting grids
      if (associated(input%properties%wannierplot%plot1d)) then
        nv = size(input%properties%wannierplot%plot1d%path%pointarray)
        if (nv < 1) then
          write(*,*)
          write(*, '("Error (wfutil_plot): Wrong plot specification!")')
          stop
        end if
        np = input%properties%wannierplot%plot1d%path%steps
        If (np < nv) then
          write(*,*)
          write(*, '("Error (wfutil_plot): Wrong plot specification!")')
          stop
        end if
  
        do ist = fst, lst
          grid( ist) = gen_1d_rgrid(input%properties%wannierplot%plot1d)
        end do
      end if
  
      if (associated(input%properties%wannierplot%plot2d)) then
        v0 = input%properties%wannierplot%plot2d%parallelogram%origin%coord
        v1 = input%properties%wannierplot%plot2d%parallelogram%pointarray(1)%point%coord
        v2 = input%properties%wannierplot%plot2d%parallelogram%pointarray(2)%point%coord
        !call r3mv( ainv, wf_centers( :, ist) + cellc - s*(v1+v2), input%properties%wannierplot%plot2d%parallelogram%origin%coord)
        !call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v1-v2), input%properties%wannierplot%plot2d%parallelogram%pointarray(1)%point%coord)
        !call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v2-v1), input%properties%wannierplot%plot2d%parallelogram%pointarray(2)%point%coord)
        do ist = fst, lst
          call r3mv( ainv, wf_centers( :, ist), v3)
          input%properties%wannierplot%plot2d%parallelogram%origin%coord              = v0 + v3 + cell - v1 - v2
          input%properties%wannierplot%plot2d%parallelogram%pointarray(1)%point%coord = v0 + v3 + cell + v1 - v2
          input%properties%wannierplot%plot2d%parallelogram%pointarray(2)%point%coord = v0 + v3 + cell - v1 + v2
          grid( ist) = gen_2d_rgrid(input%properties%wannierplot%plot2d, 0)
        end do
      end if
  
      if( associated( input%properties%wannierplot%plot3d)) then
        v0 = input%properties%wannierplot%plot3d%box%origin%coord
        v1 = input%properties%wannierplot%plot3d%box%pointarray(1)%point%coord
        v2 = input%properties%wannierplot%plot3d%box%pointarray(2)%point%coord
        v3 = input%properties%wannierplot%plot3d%box%pointarray(3)%point%coord
        s = omega**(1.d0/3.d0)
        do ist = fst, lst
          call r3mv( ainv, wf_centers( :, ist) + cellc - s*(v1+v2+v3), input%properties%wannierplot%plot3d%box%origin%coord)
          call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v1-v2-v3), input%properties%wannierplot%plot3d%box%pointarray(1)%point%coord)
          call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v2-v3-v1), input%properties%wannierplot%plot3d%box%pointarray(2)%point%coord)
          call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v3-v1-v2), input%properties%wannierplot%plot3d%box%pointarray(3)%point%coord)
          grid( ist) = gen_3d_rgrid( input%properties%wannierplot%plot3d, 0)
        end do
      end if
  
      maxnpt = 0
      do ist = fst, lst
        maxnpt = max( maxnpt, grid( ist)%npt)
      end do
      allocate( zdatatot( maxnpt, fst:lst))
      zdatatot(:,:) = zzero
      ! calculate the Wannier function on the grid
  
  
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ip, ik, evecfv, apwalm, evec, wfmt, wfir, zdata, phase, is, ia, ias)
#endif
      allocate( zdata( maxnpt))
      allocate( wfmt( lmmaxapw, nrmtmax, natmtot, nspinor))
      allocate( wfir( ngrtot, nspinor))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
      allocate( evecfv( nmatmax, nstfv))
      allocate( evec( nmatmax))
#ifdef USEOMP
!$OMP DO
#endif
      do iknr = 1, wf_kset%nkpt
        call findkptinset( vklnr( :, iknr), wf_kset, ip, ik)
#ifdef USEOMP
!$OMP CRITICAL (readevec)
#endif
        call wfhelp_getevec( ik, evecfv)
#ifdef USEOMP
!$OMP END CRITICAL (readevec)
#endif
        
        ! find the matching coefficients
        call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)
        
        call ws_weight( dble( cell), dble( cell), wf_kset%vkl( :, ik), wf_kset%ngridk, phase)
  
        do ist = fst, lst
          call zgemv( 'n', nmatmax, wf_nst, zone, &
                 evecfv( :, wf_fst:wf_lst), nmatmax, &
                 wf_transform( :, ist, ik), 1, zzero, &
                 evec, 1)
          
          ! calculate the wavefunctions for all states
          wfmt(:,:,:,:) = zzero
          wfir(:,:) = zzero
          zdata(:) = zzero
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              call wavefmt( 1, input%groundstate%lmaxapw, is, ia, wf_Gkset%ngk( 1, ik), apwalm, evec, lmmaxapw, wfmt( :, :, ias, 1))
            end do
          end do
          wfir( 1:wf_Gkset%ngk( 1, ik), 1) = evec( 1:wf_Gkset%ngk( 1, ik))/dsqrt( omega)
          call calc_zdata_rgrid( grid( ist), iknr, wfmt(:,:,:,1), wfir(:,1), zdata( 1:grid( ist)%npt), nosym=.true.)
#ifdef USEOMP
!$OMP CRITICAL (adddata)
#endif
          zdatatot( 1:grid( ist)%npt, ist) = zdatatot( 1:grid( ist)%npt, ist) + phase*zdata( 1:grid( ist)%npt)
#ifdef USEOMP
!$OMP END CRITICAL (adddata)
#endif
        end do
  
      end do        
#ifdef USEOMP
!$OMP END DO 
#endif
      deallocate( zdata, wfmt, wfir, apwalm, evecfv, evec)
#ifdef USEOMP
!$OMP END PARALLEL 
#endif
  
      zdatatot = zdatatot/wf_kset%nkpt
  
      ist = fst
      ! phase correction
      allocate( dist( maxnpt))
      do ist = fst, lst
        phi = 0.d0
        do ip = 1, grid( ist)%npt
          dist( ip) = norm2( grid( ist)%vpc( :, ip) - wf_centers( :, ist) - cellc(:))
          s = atan2( aimag( zdatatot( ip, ist)), dble( zdatatot( ip, ist)))
          !write(*,'(F23.16)') s
          if( s .gt. 0.5d0*pi) s = s - pi
          if( s .lt. -0.5d0*pi) s = s + pi
          !write(*,'(I,F13.6)') ip, s
          phi = phi + abs( zdatatot( ip, ist))*s
        end do
        phi = phi/sum( abs( zdatatot( :, ist)))
        ip = minloc( dist, 1) 
        s = atan2( aimag( zdatatot( ip, ist)), dble( zdatatot( ip, ist)))
        if( abs( phi - s) .gt. 0.5d0*pi) phi = phi + pi
        !write(*,*) ip
        !write(*,'(3F13.6)') wf_centers( :, ist) + cellc
        !write(*,'(3F13.6)') grid( ist)%vpc( :, ip)
        !write(*,'(3F13.6)') minval( grid( ist)%vpc( 1, :)), minval( grid( ist)%vpc( 2, :)), minval( grid( ist)%vpc( 3, :))
        !write(*,'(3F13.6)') maxval( grid( ist)%vpc( 1, :)), maxval( grid( ist)%vpc( 2, :)), maxval( grid( ist)%vpc( 3, :))
        phase = cmplx( cos( phi), -sin( phi), 8)
        zdatatot( :, ist) = phase*zdatatot( :, ist)
        rrange( :, ist) = (/minval( dble( zdatatot( :, ist))), maxval( dble( zdatatot( :, ist)))/)
        irange( :, ist) = (/minval( aimag( zdatatot( :, ist))), maxval( aimag( zdatatot( :, ist)))/)
      end do
  
      !----------------
      ! 1D case
      !----------------
      if (associated(input%properties%wannierplot%plot1d)) then
        ! Output
        write(*,*)
        write(*,'("Info (wfutil_plot):")')
        do ist = fst, lst
          write(fname,'("wannier1d-",i4.4,".dat")') ist
          open(77,file=trim(fname),status='Unknown',action='Write')
          do ip = 1, grid( ist)%npt
            ! path, |psi|^2, Re(psi), Im(psi) 
            write(77,'(4f16.6)') grid( ist)%vpd(ip), abs(zdatatot(ip, ist))**2, zdatatot(ip, ist)
            !write(77,'(2f16.6)') grid%vpd(ip), wkpt(ik)*nkptnr*abs(zdata(ip))**2
          end do
          close(77)
          write(*,'(" 1D Wannier function written to wannier1d-",i4.4,".dat")'), ist
          write(*,'(" real part range: ",2f13.6)') rrange( :, ist)
          write(*,'(" imag part range: ",2f13.6)') irange( :, ist)
          write(*,*)
          call delete_rgrid( grid( ist))
        end do
      end if
  
      !----------------
      ! 2D case
      !----------------
      if (associated(input%properties%wannierplot%plot2d)) then
        write(*,*)
        write(*,'("Info (wannier_plot):")')
        do ist = fst, lst
          write(fname,'("wannier2d-",i4.4,".xsf")') ist
          call str_strip(fname)
          call write_structure_xsf(fname)
          call write_2d_xsf(fname, 'module squared',   grid( ist)%boxl(1:3,:), grid( ist)%ngrid, grid( ist)%npt, abs(zdatatot( :, ist))**2)
          call write_2d_xsf(fname, 'real',             grid( ist)%boxl(1:3,:), grid( ist)%ngrid, grid( ist)%npt, dble(zdatatot( :, ist)))
          call write_2d_xsf(fname, 'imaginary',        grid( ist)%boxl(1:3,:), grid( ist)%ngrid, grid( ist)%npt, aimag(zdatatot( :, ist)))
          write(*,'(" 2D Wannier function written to wannier2d-",i4.4,".xsf")'), ist
          write(*,'(" real part range: ",2f13.6)') rrange( :, ist)
          write(*,'(" imag part range: ",2f13.6)') irange( :, ist)
          write(*,*)
          call delete_rgrid( grid( ist))
        end do
      end if
  
      !----------------
      ! 3D case
      !----------------
      if (associated(input%properties%wannierplot%plot3d)) then
        write(*,*)
        write(*,'("Info (wannier_plot):")')
        do ist = fst, lst
          write(fname,'("wannier3d-",i4.4,".xsf")') ist
          call str_strip(fname)
          call write_structure_xsf(fname)
          call write_3d_xsf(fname, 'squared modulus', grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, abs(zdatatot( :, ist))**2)
          call write_3d_xsf(fname, 'real',            grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, dble(zdatatot( :, ist)))
          call write_3d_xsf(fname, 'imaginary',       grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, aimag(zdatatot( :, ist)))
          write(*,'(" 3D Wannier function written to wannier3d-",i4.4,".xsf")'), ist
          write(*,'(" real part range: ",2f13.6)') rrange( :, ist)
          write(*,'(" imag part range: ",2f13.6)') irange( :, ist)
          write(*,*)
          !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))
  
          ! Gaussian cube-format
          write(fname,'("wannier3d-",i4.4,".cube")') ist
          call str_strip(fname)
          call write_3d_cube(fname, 'squared modulus', grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, abs(zdatatot( :, ist))**2)
          call delete_rgrid( grid( ist))
        end do
      end if
  
      deallocate( zdatatot, grid)
     
      return
    end subroutine wfutil_plot
end module mod_wannier_util
