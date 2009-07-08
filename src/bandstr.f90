! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandstr
! !INTERFACE:


subroutine bandstr
  ! !USES:
  use modinput
  use modmain
  use FoX_wxml
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
  !EOP
  !BOC
  implicit none
  ! local variables
  integer::lmax, lmmax, l, m, lm
  integer::ik, ispn, is, ia, ias, iv, ist
  real(8)::emin, emax, sum
  character(256)::fname
  ! allocatable arrays
  real(8), allocatable :: evalfv(:, :)
  real(8), allocatable :: e(:, :)
  ! low precision for band character array saves memory
  real(4), allocatable :: bc(:, :, :, :)
  complex(8), allocatable :: dmat(:, :, :, :, :)
  complex(8), allocatable :: apwalm(:, :, :, :, :)
  complex(8), allocatable :: evecfv(:, :, :)
  complex(8), allocatable :: evecsv(:, :)
  character(128)::buffer,buffer1
  type(xmlf_t), save::xf

  ! initialise universal variables
  call init0
  call init1
  ! allocate array for storing the eigenvalues
  allocate(e(nstsv, nkpt))
  ! maximum angular momentum for band character
  lmax=min(3, input%groundstate%lmaxapw)
  lmmax=(lmax+1)**2
  if (input%properties%bandstructure%character) then
     allocate(bc(0:lmax, natmtot, nstsv, nkpt))
  end if
  ! read density and potentials from file
  call readstate
  ! read Fermi energy from file
  call readfermi
  ! find the new linearisation energies
  call linengy
  ! generate the APW radial functions
  call genapwfr
  ! generate the local-orbital radial functions
  call genlofr
  ! compute the overlap radial integrals
  call olprad
  ! compute the Hamiltonian radial integrals
  call hmlrad
  emin=1.d5
  emax=-1.d5
  ! begin parallel loop over k-points
#ifdef KSMP
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(evalfv,evecfv,evecsv) &
  !$OMP PRIVATE(dmat,apwalm) &
  !$OMP PRIVATE(ispn,ist,is,ia,ias) &
  !$OMP PRIVATE(l,m,lm,sum)
  !$OMP DO
#endif
  do ik=1, nkpt
     allocate(evalfv(nstfv, nspnfv))
     allocate(evecfv(nmatmax, nstfv, nspnfv))
     allocate(evecsv(nstsv, nstsv))
     !$OMP CRITICAL
     write(*, '("Info(bandstr): ", I6, " of ", I6, " k-points")') ik, nkpt
     !$OMP END CRITICAL
     ! solve the first- and second-variational secular equations
     call seceqn(ik, evalfv, evecfv, evecsv)
     do ist=1, nstsv
        ! subtract the Fermi energy
        e(ist, ik)=evalsv(ist, ik)-efermi
        ! add scissors correction
        if (e(ist, ik).gt.0.d0) e(ist, ik)=e(ist, ik)+input%properties%bandstructure%scissor
        !$OMP CRITICAL
        emin=min(emin, e(ist, ik))
        emax=max(emax, e(ist, ik))
        !$OMP END CRITICAL
     end do
     ! compute the band characters if required
     if (input%properties%bandstructure%character) then
        allocate(dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
        allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
        ! find the matching coefficients
        do ispn=1, nspnfv
           call match(ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, ik), &
                sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
        end do
        ! average band character over spin and m for all atoms
        do is=1, nspecies
           do ia=1, natoms(is)
              ias=idxas(ia, is)
              ! generate the diagonal of the density matrix
              call gendmat(.true., .true., 0, lmax, is, ia, ngk(:, ik), apwalm, evecfv, &
                   evecsv, lmmax, dmat)
              do ist=1, nstsv
                 do l=0, lmax
                    sum=0.d0
                    do m=-l, l
                       lm=idxlm(l, m)
                       do ispn=1, nspinor
                          sum=sum+dble(dmat(lm, lm, ispn, ispn, ist))
                       end do
                    end do
                    bc(l, ias, ist, ik)=real(sum)
                 end do
              end do
           end do
        end do
        deallocate(dmat, apwalm)
     end if
     deallocate(evalfv, evecfv, evecsv)
     ! end loop over k-points
  end do
#ifdef KSMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
  emax=emax+(emax-emin)*0.5d0
  emin=emin-(emax-emin)*0.5d0
  ! output the band structure
  call xml_OpenFile ("bandstructure.xml", xf, replace=.true.,pretty_print=.true.)

  if (.not. input%properties%bandstructure%character) then
     open(50, file='BAND.OUT', action='WRITE', form='FORMATTED')
     call xml_NewElement(xf,"bandstructure")
     call xml_NewElement(xf,"title")
     call xml_AddCharacters(xf,trim(input%title))
     call xml_endElement(xf,"title")
     do ist=1, nstsv
        call xml_NewElement(xf,"band")
        do ik=1, nkpt
           write(50, '(2G18.10)') dpp1d(ik), e(ist, ik)
           call xml_NewElement(xf,"point")
           write(buffer,'(5G18.10)')dpp1d(ik)
           call xml_AddAttribute(xf, "distance", trim(adjustl(buffer)))
           write(buffer,'(5G18.10)')e(ist, ik)
           call xml_AddAttribute(xf, "eval", trim(adjustl(buffer)))
             call xml_endElement(xf,"point")
        end do
        call xml_endElement(xf,"band")
        write(50, '("     ")')
     end do
     close(50)
     write(*, *)
     write(*, '("Info(bandstr):")')
     write(*, '(" band structure plot written to BAND.OUT")')
  else
   call xml_NewElement(xf,"bandstructure")
     call xml_addAttribute(xf,"character","true")
     call xml_NewElement(xf,"title")
     call xml_AddCharacters(xf,trim(input%title))
     call xml_endElement(xf,"title")
     do is=1, nspecies
      call xml_NewElement(xf,"species")
      call xml_addAttribute(xf,"name",trim(spname(is)))
       call xml_addAttribute(xf,"chemicalSymbol",trim( input%structure%speciesarray(is)%species%chemicalSymbol))
        do ia=1, natoms(is)
          call xml_NewElement(xf,"atom")
          write(buffer,'(5G18.10)')atposc(:, ia, is)
           call xml_addAttribute(xf,"coord",trim(adjustl(buffer)))
           ias=idxas(ia, is)
           write(fname, '("BAND_S", I2.2, "_A", I4.4, ".OUT")') is, ia
           open(50, file=trim(fname), action='WRITE', form='FORMATTED')

           do ist=1, nstsv
           call xml_NewElement(xf,"band")
              do ik=1, nkpt
                 ! sum band character over l
                 sum=0.d0
                 do l=0, lmax
                    sum=sum+bc(l, ias, ist, ik)
                 end do
           call xml_NewElement(xf,"point")
           write(buffer,'(5G18.10)')dpp1d(ik)
           call xml_AddAttribute(xf, "distance", trim(adjustl(buffer)))
           write(buffer,'(5G18.10)')e(ist, ik)
           call xml_AddAttribute(xf, "eval", trim(adjustl(buffer)))
              write(buffer,'(5G18.10)')sum
           call xml_AddAttribute(xf, "sum", trim(adjustl(buffer)))
           Do l=0,lmax
           call xml_NewElement(xf,"bc")
           write(buffer,*)l
           call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
           write(buffer,'(5G18.10)')bc(l, ias, ist, ik)
           call xml_AddAttribute(xf, "character", trim(adjustl(buffer)))
           call xml_endElement(xf,"bc")
           end do
           call xml_endElement(xf,"point")
                 write(50, '(2G18.10, 8F12.6)') dpp1d(ik), e(ist, ik), sum, &
                      (bc(l, ias, ist, ik), l = 0, lmax)
              end do
               call xml_endElement(xf,"band")
              write(50, '("	  ")')
           end do
           call xml_endElement(xf,"atom")
           close(50)
        end do
        call xml_endElement(xf,"species")
     end do
     write(*, *)
     write(*, '("Info(bandstr):")')
     write(*, '(" band structure plot written to BAND_Sss_Aaaaa.OUT")')
     write(*, '("	for all species and atoms")')
  end if
  write(*, *)
  write(*, '(" Fermi energy is at zero in plot")')
  ! output the vertex location lines
  open(50, file='BANDLINES.OUT', action='WRITE', form='FORMATTED')
  do iv=1, nvp1d
 call xml_NewElement(xf,"vertex")

 write(buffer,'(5G18.10)')dvp1d(iv)
call xml_addAttribute(xf,"distance",trim(adjustl(buffer)))
 write(buffer,'(5G18.10)')emax
call xml_addAttribute(xf,"upperboundary",trim(adjustl(buffer)))
 write(buffer,'(5G18.10)')emin
call xml_addAttribute(xf,"lowerboundary",trim(adjustl(buffer)))
call xml_addAttribute(xf,"label",trim(adjustl(&
input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label)))
 write(buffer,'(5G18.10)')&
input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
call xml_addAttribute(xf,"coord",trim(adjustl(buffer)))
call xml_endElement(xf,"vertex")
     write(50, '(2G18.10)') dvp1d(iv), emin
     write(50, '(2G18.10)') dvp1d(iv), emax
     write(50, '("     ")')
  end do
  close(50)
  write(*, *)
  write(*, '(" vertex location lines written to BANDLINES.OUT")')
  write(*, *)
  deallocate(e)
  if (input%properties%bandstructure%character) deallocate(bc)
  call xml_close(xf)
  return
end subroutine bandstr
!EOC
