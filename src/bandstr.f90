! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bandstr
! !INTERFACE:
!
!
Subroutine bandstr
  ! !USES:
  use modinput
  use modmain
  use modmpi
  use FoX_wxml
  use m_write_hdf5, only: hdf5_bandstructure_output
  use precision, only: dp 


  ! !DESCRIPTION:
  !   Produces a band structure along the path in reciprocal-space which connects
  !   the vertices in the array {\tt vvlp1d}. The band structure is obtained from
  !   the second-variational eigenvalues and is written to the file {\tt BAND.OUT}
  !   with the Fermi energy set to zero. If required, band structures are plotted
  !   to files {\tt BAND\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species {\tt ss},
  !   which include the band characters for each $l$ component of that atom in
  !   columns 4 onwards. Column 3 contains the sum over $l$ of the characters.
  !   Vertex location lines are written to {\tt BANDLINES.OUT}.
  !
  ! !REVISION HISTORY:
  !   Created June 2003 (JKD)
  !   Modified June 2012 (DIN)
  !   Modified March 2014 (UW)
  !   Modified June 2018 (SeTi)
  !EOP
  !BOC
  Implicit None
  ! local variables
  Integer :: lmax, lmmax, l, m, lm
  Integer :: ik, ispn, is, ia, ias, iv, ist
  Real (8) :: emin, emax, sum
  Character (256) :: fname
  ! allocatable arrays
  Real (8), Allocatable :: evalfv (:, :)
  ! low precision for band character array saves memory
  Real (4), Allocatable :: bc (:, :, :, :)
  Complex (8), Allocatable :: dmat (:, :, :, :, :)
  Complex (8), Allocatable :: apwalm (:, :, :, :, :)
  Complex (8), Allocatable :: evecfv (:, :, :)
  Complex (8), Allocatable :: evecsv (:, :)
  Character (128) :: buffer
  Type (xmlf_t), Save :: xf

  character(:), allocatable :: label_names
  real(dp), allocatable :: lable_coordinates(:, :)

  ! initialise global variables
  Call init0
  Call init1

  !------------------------------------------------------------------
  ! In case of hybrid functionals, one is asked to use Wannier tools
  !------------------------------------------------------------------
  if (xctype(1) >= 400) then
    if (rank==0) then
      write(*,*)
      write(*,*) 'ERROR(bandstr): Please employ the WANNIER tool to interpolate the electronic band structure!'
      write(*,*)
    end if
    call terminate()
  end if

  ! maximum angular momentum for band character
  lmax = Min (3, input%groundstate%lmaxapw)
  lmmax = (lmax+1) ** 2
  If (input%properties%bandstructure%character) Then
      Allocate (bc(0:lmax, natmtot, nstsv, nkpt))
  End If
  ! read density and potentials from file
  Call readstate
  ! read Fermi energy from file
  Call readfermi
  ! find the new linearisation energies
  Call linengy
  ! generate the APW radial functions
  Call genapwfr
  ! generate the local-orbital radial functions
  Call genlofr
  ! compute the overlap radial integrals
  Call olprad
  ! compute the Hamiltonian radial integrals
  call MTNullify(mt_hscf)
  call MTInitAll(mt_hscf)
  call hmlint(mt_hscf)

  ! compute "relativistic mass"
  Call genmeffig
  emin = 1.d5
  emax = - 1.d5

  !---------------------------------------
  ! begin parallel loop over k-points
  !---------------------------------------
#ifdef MPI
  Call MPI_barrier (MPI_COMM_WORLD, ierr)
  splittfile = .True.
  Do ik = firstk(rank, nkpt), lastk(rank, nkpt)
#else
  splittfile = .False.
  Do ik = 1, nkpt
#endif
    Allocate (evalfv(nstfv, nspnfv))
    Allocate (evecfv(nmatmax, nstfv, nspnfv))
    Allocate (evecsv(nstsv, nstsv))
    ! initialise the eigenvectors if we use the Davidson eigensolver
    if (input%groundstate%solver%type.eq.'Davidson') evecfv=zzero
    ! solve the first- and second-variational secular equations
    Call seceqn (ik, evalfv, evecfv, evecsv)
    Do ist = 1, nstsv
        ! subtract the Fermi energy
        evalsv (ist, ik) = evalsv (ist, ik) - efermi
        ! add scissors correction
        If (evalsv(ist, ik) .Gt. 0.d0) evalsv (ist, ik) = evalsv (ist, ik) + &
            & input%properties%bandstructure%scissor
        emin = Min (emin, evalsv(ist, ik))
        emax = Max (emax, evalsv(ist, ik))
    End Do
    ! compute the band characters if required
    If (input%properties%bandstructure%character) Then
        Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
        Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, &
            & nspnfv))
        ! find the matching coefficients
        apwalm = zzero
        Do ispn = 1, nspnfv
          Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, &
                & ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, &
                & ispn))
        End Do
        ! average band character over spin and m for all atoms
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
              ias = idxas (ia, is)
              ! generate the diagonal of the density matrix
              Call gendmat (.True., .True., 0, lmax, is, ia, ngk(:, &
                  & ik), apwalm, evecfv, evecsv, lmmax, dmat)
              Do ist = 1, nstsv
                Do l = 0, lmax
                    sum = 0.d0
                    Do m = - l, l
                      lm = idxlm (l, m)
                      Do ispn = 1, nspinor
                          sum = sum + dble (dmat(lm, lm, ispn, &
                              & ispn, ist))
                      End Do
                    End Do
                    bc (l, ias, ist, ik) = real (sum)
                End Do
              End Do
          End Do
        End Do
        Deallocate (dmat, apwalm)
    End If ! character
    Deallocate (evalfv, evecfv, evecsv)
    ! end loop over k-points
  End Do
  call mt_hscf%release()
#ifdef MPI
  If (input%properties%bandstructure%character) Then
    Call mpi_allgatherv_ifc(nkpt, (lmax+1)*natmtot*nstsv, rlpbuf=bc)
  End If
  Call mpi_allgatherv_ifc(nkpt, nstsv, rbuf=evalsv)
  Call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif

  if (allocated(meffig)) deallocate(meffig)
  if (allocated(m2effig)) deallocate(m2effig)

  emax = emax + (emax-emin) * 0.5d0
  emin = emin - (emax-emin) * 0.5d0

  allocate(lable_coordinates(3, nvp1d))
  label_names = trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(1)%point%label))
  lable_coordinates(:, 1) = input%properties%bandstructure%plot1d%path%pointarray(1)%point%coord
  do iv=2, nvp1d
    label_names = label_names // "," // trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label))
    lable_coordinates(:, iv) = input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
  end do

  if (input%properties%bandstructure%character) then
    call hdf5_bandstructure_output(mpiglobal, 'properties.h5', '/', evalsv, [emin, emax], dpp1d, label_names, dvp1d, lable_coordinates, characters=bc)
  else
    call hdf5_bandstructure_output(mpiglobal, 'properties.h5', '/', evalsv, [emin, emax], dpp1d, label_names, dvp1d, lable_coordinates)
  end if  

  !------------------------------
  ! output the band structure
  !------------------------------
  if (rank==0) then
    Call xml_OpenFile ("bandstructure.xml", xf, replace=.True., pretty_print=.True.)
    Call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
                      &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')

    If ( .Not. input%properties%bandstructure%character) Then

      Open (50, File='BAND.OUT', Action='WRITE', Form='FORMATTED')
      Call xml_NewElement (xf, "bandstructure")
      Call xml_NewElement (xf, "title")
      Call xml_AddCharacters (xf, trim(input%title))
      Call xml_endElement (xf, "title")
      Do ist = 1, nstsv
        Call xml_NewElement (xf, "band")
        Do ik = 1, nkpt
          Write (50, '(2G18.10)') dpp1d (ik), evalsv (ist, ik)
          Call xml_NewElement (xf, "point")
          Write (buffer, '(5G18.10)') dpp1d (ik)
          Call xml_AddAttribute (xf, "distance", &
                & trim(adjustl(buffer)))
          Write (buffer, '(5G18.10)') evalsv (ist, ik)
          Call xml_AddAttribute (xf, "eval", &
                & trim(adjustl(buffer)))
          Call xml_endElement (xf, "point")
        End Do
        Call xml_endElement (xf, "band")
        Write (50, '("     ")')
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(bandstr):")')
      Write (*, '(" band structure plot written to BAND.OUT")')

    Else

      Call xml_NewElement (xf, "bandstructure")
      Call xml_AddAttribute (xf, "character", "true")
      Call xml_NewElement (xf, "title")
      Call xml_AddCharacters (xf, trim(input%title))
      Call xml_endElement (xf, "title")
      Do is = 1, nspecies
        Call xml_NewElement (xf, "species")
        Call xml_AddAttribute (xf, "name", trim(spname(is)))
        Call xml_AddAttribute (xf, "chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol))
        Do ia = 1, natoms (is)
          Call xml_NewElement (xf, "atom")
          Write (buffer, '(5G18.10)') atposc (:, ia, is)
          Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
          ias = idxas (ia, is)
          Write (fname, '("BAND_S", I2.2, "_A", I4.4, ".OUT")') is, ia
          Open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
          Do ist = 1, nstsv
            Call xml_NewElement (xf, "band")
            Do ik = 1, nkpt
              ! sum band character over l
              sum = 0.d0
              Do l = 0, lmax
                  sum = sum + bc (l, ias, ist, ik)
              End Do
              Call xml_NewElement (xf, "point")
              Write (buffer, '(5G18.10)') dpp1d (ik)
              Call xml_AddAttribute (xf, "distance", &
                    & trim(adjustl(buffer)))
              Write (buffer, '(5G18.10)') evalsv (ist, ik)
              Call xml_AddAttribute (xf, "eval", &
                    & trim(adjustl(buffer)))
              Write (buffer, '(5G18.10)') sum
              Call xml_AddAttribute (xf, "sum", &
                    & trim(adjustl(buffer)))
              Do l = 0, lmax
                  Call xml_NewElement (xf, "bc")
                  Write (buffer,*) l
                  Call xml_AddAttribute (xf, "l", &
                      & trim(adjustl(buffer)))
                  Write (buffer, '(5G18.10)') bc (l,ias,ist,ik)
                  Call xml_AddAttribute (xf, "character", trim(adjustl(buffer)))
                  Call xml_endElement (xf, "bc")
              End Do
              Call xml_endElement (xf, "point")
              Write (50, '(2G18.10, 8F12.6)') dpp1d (ik), evalsv(ist, ik), sum, (bc(l,ias,ist,ik), l=0,lmax)
            End Do
            Call xml_endElement (xf, "band")
            Write (50, '("	  ")')
          End Do
          Call xml_endElement (xf, "atom")
          Close (50)
        End Do ! ia
        Call xml_endElement (xf, "species")
      End Do ! is
      Write (*,*)
      Write (*, '("Info(bandstr):")')
      Write (*, '(" band structure plot written to BAND_Sss_Aaaaa.OUT")')
      Write (*, '("	for all species and atoms")')
    End If

    Write (*,*)
    Write (*, '(" Fermi energy is at zero in plot")')

    ! output the vertex location lines
    Open (50, File='BANDLINES.OUT', Action='WRITE', Form='FORMATTED')
    Do iv = 1, nvp1d
      Call xml_NewElement (xf, "vertex")
      Write (buffer, '(5G18.10)') dvp1d (iv)
      Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
      Write (buffer, '(5G18.10)') emax
      Call xml_AddAttribute (xf, "upperboundary", trim(adjustl(buffer)))
      Write (buffer, '(5G18.10)') emin
      Call xml_AddAttribute (xf, "lowerboundary", trim(adjustl(buffer)))
      Call xml_AddAttribute (xf, "label", trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label)))
      Write (buffer, '(5G18.10)') input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
      Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
      Call xml_endElement (xf, "vertex")
      Write (50, '(2G18.10)') dvp1d (iv), emin
      Write (50, '(2G18.10)') dvp1d (iv), emax
      Write (50, '("     ")')
    End Do
    Close (50)
    Write (*,*)
    Write (*, '(" Vertex location lines written to BANDLINES.OUT")')
    Write (*,*)
    Call xml_close (xf)
  end if ! rank

  If (input%properties%bandstructure%character) deallocate(bc)

  !---------------------------------------------------------------------------------------------------------
  ! Sorry! One more (1000+1, :-) output file for the band structure to be able to apply interpolation on it
  !---------------------------------------------------------------------------------------------------------
  if (rank==0) then
    open(50, File="bandstructure.dat", Action='Write', Form='Formatted')
    write(50,*) "# ", 1, nstsv, nkpt
    do ist = 1, nstsv
      do ik = 1, nkpt
        write(50,'(2I6, 3F12.6, 2G18.10)') ist, ik, vkl(:,ik), dpp1d(ik), evalsv(ist,ik)
      end do
    write(50,*)
    end do
    close(50)
  end if
  Return
End Subroutine bandstr
!EOC
