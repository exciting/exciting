subroutine readinp
  
  use modinput
  use modmain
  use optica
  
  implicit none
  integer    :: ik, ib

  call linmsg(166,'-','Input parameters')

! Initialization  
  call init0
  call init1

! Fermi energy
  call readfermi

!.................... Input parameters ...........................

  banmin=input%properties%nlo%banmin
  banmax=input%properties%nlo%banmax
  idel1=input%properties%nlo%delta
  dw=input%properties%nlo%dw
  emesh=input%properties%nlo%emesh
  sc=input%properties%nlo%sc
  tol=input%properties%nlo%tol
  v1=input%properties%nlo%chicomp(1)
  v2=input%properties%nlo%chicomp(2)
  v3 =input%properties%nlo%chicomp(3)
  
  write(166,*)
  write(166,*) 'min band used', banmin
  write(166,*) 'max band used', banmax
  write(166,*) 'maximal number of available states ', nstfv
  write(166,*) 'calculated the component:',v1,v2,v3 ,'of second order susceptibility'
  write(166,*) 'tolerance:', tol
  if(tol.gt.0.008) write(166,*)'ATTENTION: tolerence is too high'
  write(166,*) 'broadening:',idel1
  if (idel1.gt.0.009) then
    write(166,*) ' '
    write(166,*)'ATTENTION: broadening is quite high'
    write(166,*) ' '
  else if (idel1.gt.0.015) then
    write(166,*) ' '
    write(166,*)'ATTENTION: broadening is too high'
    write(166,*) ' '
  end if
  write(166,*) 'scissors shift:',sc

!................... fool proof ..................................

  if(tol.lt.1.d-6) then
    write(6,*) 'ERROR(nlo): Tolerence is too small!'
    stop
  end if

  if(v1.le.0.or.v2.le.0.or.v3.le.0) then
    write(6,*) 'ERROR(nlo): Components to be calculated are incorrect'
    write(6,*) 'ERROR(nlo): Components have to be more than zero'
    stop
  end if

  if(v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
    write(6,*) 'ERROR(nlo): Components to be calculated are incorrect'
    write(6,*) 'ERROR(nlo): Components have to be less than 3'
    stop
  end if
  
  if(banmax.gt.nstfv) then
    write(6,*) 'ERROR(nlo): banmax exceed maximum number of bands'
    write(6,*) 'ERROR(nlo): Check input!!'
    stop
  end if
  
  if(banmax.eq.0) then
    write(6,*) 'WARNING(nlo): banmax is not specified!'
    write(6,*) 'WARNING(nlo): Use default value equal to the maximum number of bands'
  end if
  
  if (banmin.gt.banmax) then
    write(6,*) 'ERROR(nlo): No transitions'
    write(6,*) 'ERROR(nlo): banmin is greater than banmax'
    stop
  end if

  
  allocate(noval(nkpt),nocond(nkpt),tot_ban(nkpt))
  noval(:)=0
  nocond(:)=0
  tot_ban(:)=0

  do ik = 1, nkpt
  
! ......... get the eigenvalues and occupancies from file ............
    call getevalfv(vkl(:,ik),evalfv(:,ik))

! ................... fool proof......................................

    if(evalfv(banmin,ik).gt.efermi) then
       write(6,*) 'ERROR(nlo): No transitions:'
       write(6,*) 'ERROR(nlo): Lowest band above fermi level'
       stop
    end if
    
!  ..... Counting the total number of valence and conduction bands .....
    do ib = 1, nstfv
      if(evalfv(ib,ik).le.efermi) then
        noval(ik)=noval(ik)+1
      else
        nocond(ik)=nocond(ik)+1
      end if
    end do ! ib
    tot_ban(ik)=noval(ik)+nocond(ik)
    
!............... kpt loop ends .......................................
  end do ! ik

return
end  subroutine
