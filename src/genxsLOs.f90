
subroutine genxsLOs
    use modinput
    use modmain
    use modspdeflist
    use modspdb
    use FoX_wxml
    use FoX_dom
    implicit none
    
    integer :: is, ias
    integer :: n, nn, nr, l, lmax
    real(8) :: energy, elo, ehi, flo, fhi, emi, fmi
    real(8) :: enode(0:input%groundstate%xsLO%maxnodes)
    real(8) :: vr(nrmtmax)
    real(8) :: p0s(nrmtmax), p1s(nrmtmax)
    real(8) :: q0s(nrmtmax), q1s(nrmtmax)
    real(8) :: hp0(nrmtmax)
    real(8) :: elorb(0:input%groundstate%xsLO%maxnodes, &
    &                0:input%groundstate%lmaxapw)
    
    integer :: nxslo(0:input%groundstate%lmaxapw)
    type(xmlf_t), save :: xf
    character(100) :: buffer
    integer :: io, ist, lx, ilo

    if (input%groundstate%xsLO%lmax > input%groundstate%lmaxapw) then
      input%groundstate%xsLO%lmax = input%groundstate%lmaxapw
    end if
    
    !----------------------------
    ! Generate xsLO species file
    !----------------------------

    write(*,*) '-----------------'
    write(*,*) 'Energy parameters'
    write(*,*) '-----------------'
    
    do is = 1, nspecies
    
      write(*,*) 'species',is
      lmax = max(maxval(spl(:,is)),input%groundstate%xsLO%lmax)
      
      nr = nrmt(is)
      ias = idxas(1,is)
      vr(1:nr) = veffmt(1,1:nr,ias)*y00
      elorb(:,:) = 0.d0
      do l = 0, lmax
        write(*,*) 'l=', l
        do n = 0, input%groundstate%xsLO%maxnodes
          enode(n) = 0d0
          call rdirac(0,n+l+1,l,l+1,nr,spr(:,is),vr,enode(n),p0s,q0s,.false.,.false.)
        end do ! n
        do n = 0, input%groundstate%xsLO%maxnodes
          ehi = enode(n)
          call rschroddme(0,l,0,ehi,nr,spr(:,is),vr,nn,p0s,hp0,q0s,q1s)
          fhi = hp0(nrmt(is))
          if (p0s(nrmt(is))==p0s(nrmt(is)-1)) then
            elo=ehi
          else
            if (n==0) then 
              elo = 2*enode(0)-enode(1) ! assuming lowest eigenenergy is  negative
            else
              elo = enode(n-1)
            end if
            call rschroddme(0,l,0,elo,nr,spr(:,is),vr,nn,p0s,hp0,q0s,q1s) 
            flo = hp0(nrmt(is))
            if (ehi<elo) then
              write(*,*) 'Oops! Error(genxsLOs): ehi<elo'
              stop
            end if
            do while (ehi-elo>1d-6)
              emi = 0.5d0*(ehi+elo)
              call rschroddme(0,l,0,emi,nr,spr(:,is),vr,nn,p0s,hp0,q0s,q1s)
              fmi = hp0(nrmt (is))
              if (fmi*fhi<0) then
                flo = fmi
                elo = emi
              else
                fhi = fmi
                ehi = emi
              end if
            end do
          end if
	        elorb(n,l) = 0.5d0*(enode(n)+0.5d0*(ehi+elo))
          write(*,*) 'n=', n, elorb(n,l)
        end do ! n
        write(*,*)
      end do ! l
  
      call xml_OpenFile (trim(spsymb(is))//'-xsLO.xml', xf, replace=.true.,pretty_print=.true.)
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
      ! DOUBLE radial grid points
      write(buffer,*) 2*nrmt(is)
      call xml_AddAttribute(xf, "radialmeshPoints", trim(adjustl(buffer)))
      call xml_endElement(xf, "muffinTin")

      do ist = 1, spnst(is)
        call xml_NewElement(xf,"atomicState")
        write(buffer,*) spn(ist,is) 
        call xml_AddAttribute(xf,"n",trim(adjustl(buffer)))
        write(buffer,*) spl(ist,is)
        call xml_AddAttribute(xf,"l",trim(adjustl(buffer)))
        write(buffer,*) spk(ist, is)
        call xml_AddAttribute(xf,"kappa",trim(adjustl(buffer)))
        write(buffer,'(G14.6)') spocc(ist,is)
        call xml_AddAttribute(xf,"occ",trim(adjustl(buffer)))
        if (spcore(ist,is)) then
          buffer="true"
        else
          buffer="false"
        endif
        call xml_AddAttribute(xf, "core", trim(adjustl(buffer)))
        call xml_endElement(xf, "atomicState")
      end do

      call xml_NewElement(xf,"basis")
      
      ! Default augmentation
      call xml_NewElement(xf,"default")
      call xml_AddAttribute(xf,"type","lapw")
      if (size(speziesdeflist(is)%sp%basis%default%wfarray)>0) then
        write(buffer,'(F8.4)') speziesdeflist(is)%sp%basis%default%wfarray(1)%wf%trialEnergy
      else
        write(buffer,'(F8.4)') speziesdeflist(is)%sp%basis%default%trialEnergy
      end if
      call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
      call xml_AddAttribute(xf,"searchE","false")
      call xml_endElement(xf,"default")
      
      nxslo(:) = 0
      do ist = 1, spnst(is)
        l = spl(ist,is)
	      nxslo(l) = max(spn(ist,is)-l,nxslo(l))
      end do
      
      do l = 0, lmax

        ! Augmentation type      
        call xml_AddNewline(xf)
        call xml_NewElement(xf,"custom")
	      write(buffer,*) l
        call xml_AddAttribute(xf,"l",trim(adjustl(buffer)))
        
if (.false.) then        
        if (apword(l,is)==1) then
          call xml_AddAttribute(xf,"type","apw")
          write(buffer,'(F8.4)') apwe0(1,l,is)
          call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
          call xml_AddAttribute(xf,"searchE","false")  
        else if (apword(l,is)==2) then
          call xml_AddAttribute(xf,"type","lapw")
          write(buffer,'(F8.4)') apwe0(1,l,is)
          call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
          call xml_AddAttribute(xf,"searchE","false") 
        else
          do io = 1, apword(l,is)
            call xml_NewElement(xf,"wf")
            write(buffer,*) apwdm(io,l,is)
            call xml_AddAttribute(xf,"matchingOrder",trim(adjustl(buffer)))
            write(buffer,'(F8.4)') apwe0(io,l,is)
            call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", "false")
            call xml_endElement(xf, "wf")
          end do
        end if
else
        ! LAPW for all
        call xml_AddAttribute(xf,"type","lapw")
        write(buffer,'(F8.4)') apwe0(1,l,is)
        call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
        call xml_AddAttribute(xf,"searchE","false")
end if
        call xml_endElement(xf,"custom")
      
        ! valence electrons LO's
        do ilo = 1, nlorb(is)
          lx = lorbl(ilo,is)
          if (lx==l) then
            call xml_NewElement(xf,"lo")
            write(buffer,*) l
            call xml_AddAttribute(xf,"l",trim(adjustl(buffer)))
            do io = 1, lorbord(ilo,is)
              call xml_NewElement(xf,"wf")
              write(buffer,*) lorbdm(io,ilo,is)
              ! ED-lo is present -> skip xsLO with similar number of nodes
              if (lorbdm(io,ilo,is)>1) nxslo(l) = nxslo(l)+lorbdm(io,ilo,is)-1
              ! ----
              call xml_AddAttribute(xf,"matchingOrder",trim(adjustl(buffer)))
              write(buffer,'(F8.4)') lorbe0(io,ilo,is)
              call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
              if (lorbve(io,ilo,is)) then
                call xml_AddAttribute(xf,"searchE","true")
              else
                call xml_AddAttribute(xf,"searchE","false")
              end if
              call xml_endElement(xf,"wf")
            end do ! io
            call xml_endElement(xf,"lo")
          end if
        end do ! ilo
              
        !------ lapw+lo for all excited states
        if (apword(l,is)==2) then
          call xml_NewElement(xf,"lo")                                      
          write(buffer,*) l                                                 
          call xml_AddAttribute(xf,"l",trim(adjustl(buffer)))               
          do io = 0, 1
            call xml_NewElement(xf,"wf")                                    
            write(buffer,*) io                                              
            call xml_AddAttribute(xf,"matchingOrder",trim(adjustl(buffer)))
            write(buffer,'(F8.4)') apwe0(1,l,is)
            call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
            call xml_AddAttribute(xf,"searchE","false")
            call xml_endElement(xf,"wf")
          end do
          call xml_endElement(xf,"lo")
          ! therefore one needs to skip the first xsLO
          nxslo(l) = nxslo(l)+1
        end if

        ! xs LO's
        if (l<=input%groundstate%xsLO%lmax) then
          do n = nxslo(l), input%groundstate%xsLO%maxnodes
            if (elorb(n,l) <= input%groundstate%xsLO%emax) then
              call xml_NewElement(xf,"lo")                                      
              write(buffer,*) l                                                 
              call xml_AddAttribute(xf,"l",trim(adjustl(buffer)))               
              
              !------ sc-LO style
              io = 0
              call xml_NewElement(xf,"wf")                                    
              write(buffer,*) io                                              
              call xml_AddAttribute(xf,"matchingOrder",trim(adjustl(buffer)))
              write(buffer,'(F8.4)') apwe0(1,l,is)
              call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
              call xml_AddAttribute(xf,"searchE","false")
              call xml_endElement(xf,"wf")
              
              call xml_NewElement(xf,"wf")
              write(buffer,*) io                                              
              call xml_AddAttribute(xf,"matchingOrder",trim(adjustl(buffer))) 
              write(buffer,'(F8.4)') elorb(n,l)                               
              call xml_AddAttribute(xf,"trialEnergy",trim(adjustl(buffer)))
              call xml_AddAttribute(xf,"searchE","false")                     
              call xml_endElement(xf,"wf")                                    
              
              call xml_endElement(xf,"lo")
            end if                                                              
          end do ! n                                                            
        end if

      end do ! l
      
      call xml_endElement(xf, "basis")

      call xml_AddNewline(xf)
      call xml_AddNewline(xf)
      call xml_AddComment(xf," This file was automatically generated ")
      call xml_AddNewline(xf)
      call xml_AddNewline(xf)

      call xml_endElement (xf, "sp")
      call xml_endElement (xf, "spdb")
      call xml_Close(xf)

    end do ! is  
    
    return    
end subroutine
