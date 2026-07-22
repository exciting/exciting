module species_file_update

  use mod_atoms, only: nspecies, spsymb, spname, spzn, spmass, sprmin, sprmax, spnst, spn, spl, spk, spocc, spcore, idxas
  use mod_muffin_tin, only: rmt, nrmt
  use mod_APW_LO, only: maxlapw, apword, apwe, apwve, nlorb, lorbl, lorbord, lorbdm, lorbe, lorbve, apwdm
  use modinput, only: input
  use modspdeflist
  use precision, only: i32
  use FoX_wxml
  use FoX_dom

  implicit none

  private 
  public :: update_species

contains

  !> This subroutine writes a new species file with the extension _scf after every exciting calculation.
  !> It conatains the same basis functions as the initial species file, but linearization energies are
  !> replaced by the lineariztion energies used in the last scf iteration.
  subroutine update_species()

    ! local variables
    type(xmlf_t), save :: xf
    character(100) :: buffer
    integer(i32) :: is, io, ist, lx, ilo, iapw, ias

    !---------------------------------------------------------
    ! Write down new (updated) xml species file
    !---------------------------------------------------------
    do is = 1, nspecies
       ias = idxas(1, is)

       call xml_OpenFile (trim(spsymb(is))//'_scf.xml', xf, replace=.true.,pretty_print=.true.)
       write(buffer,*) trim(input%groundstate%findlinentype)
       call xml_AddComment(xf, &
            " This file was automatically generated."//char(10)// &
            " Linearization energies with searchE=true have been replaced by energies"//char(10)// &
            " found using findlinentype='"//trim(adjustl(buffer))//"' method."//char(10)// &
            " If a principal quantum number was given, linearization energies have been"//char(10)// &
            " set to energies found using the Wigner-Seitz algorithm in the first scf iteration.")
       call xml_NewElement (xf, "spdb")
       call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XMLSchema-instance", "xsi")
       call xml_AddAttribute(xf, "xsi:noNamespaceSchemaLocation", "../../xml/species.xsd" )
       call xml_NewElement (xf, "sp")
       call xml_AddAttribute(xf, "chemicalSymbol", trim(adjustl(spsymb(is))) )
       call xml_AddAttribute(xf, "name", trim(adjustl(spname(is))) )
       write(buffer,'(G14.6)') spzn(is)
       call xml_AddAttribute(xf, "z", trim(adjustl(buffer)))
       write(buffer,'(G18.10)') spmass(is)
       call xml_AddAttribute(xf, "mass", trim(adjustl(buffer)))
       call xml_NewElement (xf, "muffinTin")
       write(buffer,'(G14.6)') sprmin(is)
       call xml_AddAttribute(xf, "rmin", trim(adjustl(buffer)))
       write(buffer,'(F10.4)') rmt(is)
       call xml_AddAttribute(xf, "radius", trim(adjustl(buffer)))
       write(buffer,'(F10.4)') sprmax(is)
       call xml_AddAttribute(xf, "rinf", trim(adjustl(buffer)))
       write(buffer,*) nrmt(is)
       call xml_AddAttribute(xf, "radialmeshPoints", trim(adjustl(buffer)))
       call xml_endElement(xf, "muffinTin")

       do ist = 1, spnst(is)
          call xml_NewElement (xf, "atomicState")
          write(buffer,*) spn(ist, is) 
          call xml_AddAttribute(xf, "n", trim(adjustl(buffer)))
          write(buffer,*) spl(ist, is)
          call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
          write(buffer,*) spk(ist, is)
          call xml_AddAttribute(xf, "kappa", trim(adjustl(buffer)))
          write(buffer,'(G14.6)') spocc(ist, is)
          call xml_AddAttribute(xf, "occ", trim(adjustl(buffer)))
          if (spcore(ist,is)) then
             buffer="true"
          else
             buffer="false"
          endif
          call xml_AddAttribute(xf, "core", trim(adjustl(buffer)))
          call xml_endElement(xf, "atomicState")
       end do !ist

       !--------------------------------------------------------------
       ! BASIS
       !--------------------------------------------------------------
       call xml_NewElement (xf, "basis")

       ! Default is always required, however all APWs up to lmaxapw are
       ! over written by the custom elements.
       call xml_NewElement (xf, "default")
       buffer=trim(speziesdeflist(is)%sp%basis%default%type)
       call xml_AddAttribute(xf, "type", trim(adjustl(buffer)))
       write(buffer,'(F8.4)') speziesdeflist(is)%sp%basis%default%trialEnergy
       call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
       if (speziesdeflist(is)%sp%basis%default%searchE) then
          call xml_AddAttribute(xf, "searchE", "true")
       else
          call xml_AddAttribute(xf, "searchE", "false")
       end if
       call xml_endElement (xf, "default")


       !   Custom augmentation type
       do iapw = 1, maxlapw
          lx = iapw - 1
          if (apword(lx, is) /= 0) then
             call xml_NewElement (xf, "custom")
             if (apword(lx, is) == 1) then
                call xml_AddAttribute(xf, "type", "apw")
             else if (apword(lx, is) == 2) then
                call xml_AddAttribute(xf, "type", "lapw")
             end if
             write(buffer,*) lx
             call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
             write(buffer,'(F8.4)') apwe(1, lx, ias)
             call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
             call xml_AddAttribute(xf, "searchE", "false")
             call xml_endElement (xf, "custom")
          endif
       end do ! iapw

       !   Local Orbitals
       do ilo = 1, nlorb(is)
          call xml_NewElement (xf, "lo")
          write(buffer,*) lorbl(ilo, is)
          call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
          do io = 1, lorbord(ilo, is)
             call xml_NewElement (xf, "wf")
             write(buffer,*) lorbdm(io, ilo, is)
             call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
             write(buffer,'(F8.4)') lorbe(io, ilo, ias)
             call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
             call xml_AddAttribute(xf, "searchE", "false")
             call xml_endElement (xf, "wf")
          end do ! io
          call xml_endElement (xf, "lo")
       end do ! ilo

       call xml_endElement (xf, "basis")
       call xml_endElement (xf, "sp")
       call xml_endElement (xf, "spdb")

       call xml_Close(xf)

    end do ! is

  end subroutine update_species
end module species_file_update
