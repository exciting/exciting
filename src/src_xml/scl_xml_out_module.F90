#include "../version.inc"
module scl_xml_out_Module
  use FoX_dom
  use mod_energy
  use mod_convergence
  use mod_eigenvalue_occupancy
  use mod_charge_and_moment
  use mod_atoms
  use mod_timing
  use mod_force
  use mod_spin
  use modmpi ,only:rank

  use modinput
  implicit none
  type(Node),  pointer :: sclDoc, root, np,npatt,energies,niter,nnewline,charges,atom,xst,timing,nscl,ngroundstate
  type(DOMConfiguration),pointer :: configo
  real(8)::scltime0=0
  character(512)::buffer
  character::newline
contains

  subroutine scl_xml_out_create
    character(10)::dat, tim
    if(rank.eq.0) then
       ! Create a new document and get a pointer to the root element, this gives you the minimum empty dom
       sclDoc => createDocument(getImplementation(), "", "info", null())
       configo => getDomConfig(scldoc)
       root => getDocumentElement(sclDoc)
       xst=>createProcessingInstruction(scldoc, "xml-stylesheet",&
            'href="'//trim(input%xsltpath)//'/info.xsl" type="text/xsl"')
       dummy => insertBefore(scldoc, xst, root)
       newline=ACHAR(10)
       nnewline =>createTextNode(scldoc, newline//" " )
       dummy => appendChild(root, nnewline)
       ngroundstate => createElementNS(sclDoc, "", "groundstate")
       dummy => appendChild(root, ngroundstate)
       nscl => createElementNS(sclDoc, "", "scl")
       dummy => appendChild(ngroundstate, nscl)
       nnewline =>createTextNode(scldoc, newline//" " )
       dummy => appendChild(root, nnewline)
       call date_and_time(date=dat, time=tim)

       write(buffer, '(A4, "-", A2, "-", A2)') dat(1:4), dat(5:6), &
            dat(7:8)
       call setAttribute(root, "date", trim(adjustl(buffer)))
       write(buffer, '( A2, ":", A2, ":", A2)') tim(1:2), tim(3:4), &
            tim(5:6)
       call setAttribute(root, "time", trim(adjustl(buffer)))
       buffer=GITHASH
        call setAttribute(root, "versionhash", trim(adjustl(buffer)))
		call setAttribute(root, "title", trim(adjustl(input%title)))

    endif

  end subroutine scl_xml_out_create
  subroutine scl_iter_xmlout()

    implicit none
    integer::is,ia,ias
    real(8)::scltime
    if(rank.eq.0) then
       nnewline =>createTextNode(scldoc,  newline//"  ")
       dummy => appendChild(nscl, nnewline)
       niter => createElementNS(sclDoc, "", "iter")
       dummy => appendChild(nscl, niter)
       write(buffer,*)iscl
       call setAttribute(niter, "iteration", trim(adjustl(buffer)))
       write(buffer,*)currentconvergence
       call setAttribute(niter, "rms", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')log10(currentconvergence)
       call setAttribute(niter, "rmslog10", trim(adjustl(buffer)))
       write(buffer,'(G22.12)') fermidos
       call setAttribute(niter, "fermidos", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')efermi
       nnewline =>createTextNode(scldoc,  newline//"   ")
       dummy => appendChild(niter, nnewline)
       energies => createElementNS(sclDoc, "", "energies")
       dummy => appendChild(niter, energies)
       write(buffer,*)engytot
       call setAttribute(energies, "totalEnergy", trim(adjustl(buffer)))
       write(buffer,*)efermi
       call setAttribute(energies, "fermiEnergy", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')evalsum
       call setAttribute(energies, "sum-of-eigenvalues", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engykn
       call setAttribute(energies, "electronic-kinetic", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engykncr
       call setAttribute(energies, "core-electron-kinetic", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engycl
       call setAttribute(energies, "Coulomb", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engyvcl
       call setAttribute(energies, "Coulomb-potential", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engynn
       call setAttribute(energies, "nuclear-nuclear", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engyen
       call setAttribute(energies, "electron-nuclear", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engyhar
       call setAttribute(energies, "Hartree", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')engymad
       call setAttribute(energies, "Madelung", trim(adjustl(buffer)))
       nnewline =>createTextNode(scldoc,  newline//"   ")
       dummy => appendChild(niter, nnewline)
       charges => createElementNS(sclDoc, "", "charges")
       dummy => appendChild(niter, charges)
       write(buffer,'(G22.12)')chgcr
       call setAttribute(charges, "core", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')chgcrlk
       call setAttribute(charges, "core_leakage", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')chgval
       call setAttribute(charges, "valence", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')chgir
       call setAttribute(charges, "interstitial", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')chgcalc
       call setAttribute(charges, "totalcharge", trim(adjustl(buffer)))
       write(buffer,'(G22.12)') chgmttot
       call setAttribute(charges, "muffin-tin-total", trim(adjustl(buffer)))
       if (input%groundstate%chgexs.ne.0.d0) then
          write(buffer,'(G22.12)') input%groundstate%chgexs
          call setAttribute(charges, "excess", trim(adjustl(buffer)))
       end if
       do is=1, nspecies
          do ia=1, natoms(is)
             nnewline =>createTextNode(scldoc,  newline//"    ")
	     dummy => appendChild(charges, nnewline)
             atom => createElementNS(sclDoc, "", "atom")
	     dummy => appendChild(charges, atom)
             ias=idxas(ia, is)
             write(buffer,*)spsymb(is)
             call setAttribute(atom, "species", trim(adjustl(buffer)) )
             write(buffer,'(G22.12)')chgmt(ias)
             call setAttribute(atom, "muffin-tin", trim(adjustl(buffer)))
          end do
       end do

       nnewline =>createTextNode(scldoc,  newline//"   ")
       dummy => appendChild(niter, nnewline)
       timing => createElementNS(sclDoc, "", "timing")
       dummy => appendChild(niter, timing)
       call timesec(scltime)
       write(buffer,'(G22.12)')scltime-scltime0
       scltime0= scltime
       if(iscl .ge. 2)call setAttribute(timing, "itertime", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timeinit+timemat+timefv+timesv+timerho+timepot+timefor
       call setAttribute(timing, "timetot", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timeinit
       call setAttribute(timing, "timeinit", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timemat
       call setAttribute(timing, "timemat", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timefv
       call setAttribute(timing, "timefv", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timerho
       call setAttribute(timing, "timerho", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timepot
       call setAttribute(timing, "timepot", trim(adjustl(buffer)))
       write(buffer,'(G22.12)')timefor
       call setAttribute(timing, "timefor", trim(adjustl(buffer)))
       nnewline =>createTextNode(scldoc,  newline//"    ")
       dummy => appendChild(charges, nnewline)
       nnewline =>createTextNode(scldoc,  newline//"   ")
       dummy => appendChild(niter, nnewline)

    endif
  end subroutine scl_iter_xmlout
  subroutine structure_xmlout
    use mod_lattice
    use mod_constants
    implicit none
    integer::is,ia,ias,i
    type(Node),  pointer :: structure,crystal,nbasevect,nreziprvect,species,atom,forces,text,force
    if(rank.eq.0) then
       structure => createElementNS(sclDoc, "", "structure")
       dummy => appendChild(nscl, structure)
       crystal => createElementNS(sclDoc, "", "crystal")
       dummy => appendChild(structure, crystal)
       write(buffer,'(G18.10)')omega
       call setAttribute(crystal, "unitCellVolume", trim(adjustl(buffer)))
       write(buffer,'(G18.10)')(twopi**3)/omega
      call setAttribute(crystal, "BrillouinZoneVolume", trim(adjustl(buffer)))
       do i=1,3
          nbasevect => createElementNS(sclDoc, "", "basevect")
          dummy => appendChild( crystal,nbasevect)
          write(buffer,'(3G18.10)')input%structure%crystal%basevect(:,i)
          text=>createTextNode(scldoc, trim(adjustl(buffer)))
          dummy => appendChild(nbasevect,text)
       end do
         nnewline =>createTextNode(scldoc,  newline//"     ")
             dummy => appendChild(crystal, nnewline)
        do i=1,3
           nreziprvect => createElementNS(sclDoc, "", "reziprvect")
           write(buffer, '(3G18.10)') bvec(:, i)
           text=>createTextNode(scldoc, trim(adjustl(buffer)) )
           dummy => appendChild(nreziprvect, text)
           dummy => appendChild(crystal, nreziprvect)
        end do
       do is=1,nspecies
          species => createElementNS(sclDoc, "", "species")
          dummy => appendChild(structure, species)
	      nnewline =>createTextNode(scldoc,  newline//"    ")
          dummy => appendChild(structure, nnewline)
	  write(buffer,*)trim(input%structure%speciesarray(is)%species%chemicalSymbol)
	  call setAttribute(species, "chemicalSymbol", trim(adjustl(buffer)))
          do ia=1,natoms(is)
             atom => createElementNS(sclDoc, "", "atom")
             dummy => appendChild(species, atom)
             nnewline =>createTextNode(scldoc,  newline//"     ")
             dummy => appendChild(species, nnewline)
             call setcoord(atom,input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord)
             if (input%groundstate%tforce) then
                forces=> createElementNS(sclDoc, "", "forces")
                dummy => appendChild( atom,forces)
                ias=idxas(ia, is)
                force=> createElementNS(sclDoc, "", "Hellmann-Feynman")
                dummy => appendChild(forces,force)
                call setcoord(force,forcehf(:, ias))
                force=> createElementNS(sclDoc, "", "core-correction")
                dummy => appendChild(forces,force)
                call setcoord(force,forcecr(:, ias))
                force=> createElementNS(sclDoc, "", "IBS")
                dummy => appendChild(forces,force)
                call setcoord(force,forceibs(:, ias))
                force=> createElementNS(sclDoc, "", "totalforce")
                dummy => appendChild(forces,force)
                call setcoord(force,forcetot(:, ias))
                write(buffer,'(G22.12)')sqrt(forcetot(1, ias)**2+forcetot(2, ias)**2+forcetot(3, ias)**2)
                call setAttribute(forces, "Magnitude", trim(adjustl(buffer)))
                nnewline =>createTextNode(scldoc,  newline//"       ")
                dummy => appendChild(atom, nnewline)
                nnewline =>createTextNode(scldoc,  newline//"        ")
                dummy => appendChild(forces, nnewline)
             endif
          end do
          if (input%groundstate%tforce) then
             write(buffer,'(G22.12)')forcemax
             call setAttribute(structure, "forceMax", trim(adjustl(buffer)))
          endif
          nnewline =>createTextNode(scldoc,  newline//"      ")
          dummy => appendChild(species, nnewline)
       end do
       nnewline =>createTextNode(scldoc,  newline//"    ")
       dummy => appendChild(structure, nnewline)
    endif
  end subroutine structure_xmlout
  subroutine setcoord(elementnode,coord)
    type(node),pointer,intent(in)::elementnode
    real(8),intent(in)::coord(3)
    write(buffer,'(G22.12)')coord(1)
    call setAttribute(elementnode, "x", trim(adjustl(buffer)))
    write(buffer,'(G22.12)')coord(2)
    call setAttribute(elementnode, "y", trim(adjustl(buffer)))
    write(buffer,'(G22.12)')coord(3)
    call setAttribute(elementnode, "z", trim(adjustl(buffer)))
  end subroutine setcoord
  subroutine setcoorddim(elementnode,coord,dim)
    implicit none
    type(node),pointer,intent(in)::elementnode
    real(8),intent(in)::coord(3)
    integer,intent(in)::dim
    if(dim.gt.0) then
       write(buffer,'(G22.12)')coord(1)
       call setAttribute(elementnode, "x", trim(adjustl(buffer)))
       if(dim.gt.1) then
          write(buffer,'(G22.12)')coord(2)
          call setAttribute(elementnode, "y", trim(adjustl(buffer)))
          if(dim.gt.2) then
             write(buffer,'(G22.12)')coord(3)
             call setAttribute(elementnode, "z", trim(adjustl(buffer)))
          endif
       endif
    endif
  end subroutine setcoorddim
  subroutine scl_xml_write_moments()

    type(Node),pointer:: moments,moment
    integer::is,ia,ias
    if(rank.eq.0) then
       moments => createElementNS(sclDoc, "", "moments")
       dummy => appendChild(niter, moments)
       nnewline =>createTextNode(scldoc,  newline//"     ")
       dummy => appendChild(moments, nnewline)
       moment=>createElementNS(sclDoc, "", "interstitial")
       dummy => appendChild(moments, moment)
       call setcoorddim(moment,momir(1:ndmag),ndmag)
       nnewline =>createTextNode(scldoc,  newline//"     ")
       dummy => appendChild(moments, nnewline)

       moment=>createElementNS(sclDoc, "", "mommttot")
       dummy => appendChild(moments, moment)
       call setcoorddim(moment,mommttot(1:ndmag),ndmag)
       nnewline =>createTextNode(scldoc,  newline//"     ")
       dummy => appendChild(moments, nnewline)

       moment=>createElementNS(sclDoc, "", "momtot")
       dummy => appendChild(moments, moment)
       call setcoorddim(moment,momtot(1:ndmag),ndmag)
       nnewline =>createTextNode(scldoc,  newline//"     ")
       dummy => appendChild(moments, nnewline)
       do is=1, nspecies
          do ia=1, natoms(is)
             ias=idxas(ia, is)
             atom  => createElementNS(sclDoc, "", "atom")
	     dummy => appendChild(moments, atom)
	     moment=> createElementNS(sclDoc, "", "mommt")
             dummy => appendChild(atom, moment)
	     call setcoorddim(moment,mommt(1:ndmag,ias),ndmag)
             write(buffer,*)spsymb(is)
             call setAttribute(atom, "species", trim(adjustl(buffer)))
          end do
          nnewline =>createTextNode(scldoc,  newline//"    ")
          dummy => appendChild(moments, nnewline)
       end do
       nnewline =>createTextNode(scldoc,  newline//"  ")
       dummy => appendChild(niter, nnewline)
       nnewline =>createTextNode(scldoc,  newline//" ")
       dummy => appendChild(nscl, nnewline)
    endif
  end subroutine scl_xml_write_moments
  subroutine scl_xml_out_close()!
    if(rank.eq.0) then
       call destroy(sclDoc)
    endif
  end subroutine scl_xml_out_close

  subroutine scl_xml_setGndstateStatus(status)

    character(len=*)::status
    if(rank.eq.0) then
       call setAttribute(ngroundstate, "status", status)
    endif
  end subroutine scl_xml_setGndstateStatus

  subroutine scl_xml_out_write()
    if(rank.eq.0) then
       call serialize(sclDoc,"info.xml")
    endif
  end subroutine scl_xml_out_write

end module scl_xml_out_Module
