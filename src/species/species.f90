
program species

  use modspecies
  use mod_muffin_tin, only: epsedirac, epspotatom
  implicit none
  
  integer :: nn
  real(8), allocatable :: p0(:), p1(:), q0(:), q1(:)
  real(8) :: e, t0, t, t1, t2, cutoff
  character(80) :: radialgridtype

!===============================================================================
!
! Parameter set corresponding to original source code for species-program.
! Note that for some of the shipped species files the cutoff parameters are
! different (epsedirac=1.d-10 and epspotatom=eps=1.d-7)
  epsedirac=1.d-11
  epspotatom=1.d-6
  ! set default values
  inspecies=''
  ecvcut=-3.5d0
  esccut=-0.5d0
  apwdescr=''
  suffix=''
  apword=1
  apwdm(:)=0
  apwve(:)=.false.
  apwordx=0
  apwdmx(:)=0
  apwve(:)=.false.
  locorb=.true.
  locorbsc=.true.
  searchlocorb=.false.
  fullsearchlocorbsc=.false.
  
  locorbxs=.false.
  searchlocorbxs=.false.
  exsmax=30.d0
  exscut=0.35d0
  
!===============================================================================

  ! get generation strategy from input file 'species.input'
  open(50,file='species.input',action='READ',status='OLD',form='FORMATTED', &
   iostat=iostat)
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(species): error opening species.input")')
    write(*,*)
    stop
  end if
  100 continue
  read(50,*,end=300) bname
  ! check for a comment
  if ((scan(trim(bname),'!').eq.1).or.(scan(trim(bname),'#').eq.1)) goto 100
  select case(trim(bname))
  case('epsedirac')
    read(50,*,err=200) epsedirac
    if (epsedirac.lt.0.d0) then
      write(*,*)
      write(*,'("Error(species): epsedirac < 0 : ",g18.10)') epsedirac
      write(*,*)
      stop
    end if
  case('epspotatom')
    read(50,*,err=200) epspotatom
    if (epspotatom.lt.0.d0) then
      write(*,*)
      write(*,'("Error(species): epspotatom < 0 : ",g18.10)') epspotatom
      write(*,*)
      stop
    end if
  case('inspecies')
    read(50,*,err=200) inspecies
  case('ecvcut')
    read(50,*,err=200) ecvcut
  case('esccut')
    read(50,*,err=200) esccut
  case('apwdescr')
    read(50,*,err=200) apwdescr
  case('suffix')
    read(50,*,err=200) suffix
  case('apw')
    read(50,*,err=200) apword
    if ((apword.lt.1).or.(apword.gt.maxapword)) then
      write(*,*)
      write(*,'("Error(species): apword < 1 : ",A)') apword
      write(*,*)
      stop
    end if
    if (apword.gt.maxapword) then
      write(*,*)
      write(*,'("Error(species): apword too large : ",I8)') apword
      write(*,*)
      stop
    end if
    do io=1,apword
      read(50,*) apwdm(io),apwve(io)
      if (apwdm(io).lt.0) then
        write(*,*)
        write(*,'("Error(species): apwdm < 0 : ",I8)') apwdm(io)
        write(*,'(" for order ",I4)') io
        write(*,*)
        stop
      end if
    end do
  case('apwx')
    read(50,*,err=200) apwordx
    if ((apwordx.lt.1).or.(apwordx.gt.maxapword)) then
      write(*,*)
      write(*,'("Error(species): apwordx < 1 : ",A)') apwordx
      write(*,*)
      stop
    end if
    if (apwordx.gt.maxapword) then
      write(*,*)
      write(*,'("Error(species): apwordx too large : ",I8)') apwordx
      write(*,*)
      stop
    end if
    do io=1,apwordx
      read(50,*) apwdmx(io),apwvex(io)
      if (apwdmx(io).lt.0) then
        write(*,*)
        write(*,'("Error(species): apwdmx < 0 : ",I8)') apwdmx(io)
        write(*,'(" for order ",I4)') io
        write(*,*)
        stop
      end if
    end do
  case('locorb')
    read(50,*,err=200) locorb
  case('locorbsc')
    read(50,*,err=200) locorbsc
  case('searchlocorb')
    read(50,*,err=200) searchlocorb
  case('fullsearchlocorbsc')
    read(50,*,err=200) fullsearchlocorbsc
  case('')
    goto 100
  !-------------------------------------------------------------------  
  case('locorbxs')
    read(50,*,err=200) locorbxs
  case('searchlocorbxs')
    read(50,*,err=200) searchlocorbxs
  case('lxsmax')
    read(50,*,err=200) lxsmax
  case('exsmax')
    read(50,*,err=200) exsmax
  case('exscut')
    read(50,*,err=200) exscut
  !-------------------------------------------------------------------  
  case default
    write(*,*)
    write(*,'("Error(species): invalid block name : ",A)') trim(bname)
    write(*,*)
    stop
  end select
  goto 100
  200 continue
  write(*,*)
  write(*,'("Error(species): error reading from species.in")')
  write(*,'("Problem occurred in ''",A,"'' block")') trim(bname)
  write(*,*)
  stop
  300 continue
  close(50)
  
  !===============================================================================

  open(40,file='species.dat',action='READ',status='OLD',form='FORMATTED')
  10 continue
  read(40,*,iostat=iostat) nz
  if (iostat.ne.0) stop
  read(40,*) spsymb
  read(40,*) spname
  read(40,*) spmass
  read(40,*) rmt
  read(40,*) spnst
  if (spnst.gt.maxspst) then
    write(*,*)
    write(*,'("Error(species): too many states for species ",A)') trim(spname)
    write(*,*)
    stop
  end if
  do ist=1,spnst
    read(40,*) spn(ist),spl(ist),spk(ist),i
    if (ist.ge.2) then
      if (spn(ist).lt.spn(ist-1)) then
        write(*,*)
        write(*,'("Error(species): states improperly ordered")')
        write(*,'(" for species ",A)') trim(spname)
        write(*,*)
        stop
      end if
    end if
    spocc(ist)=dble(i)
  end do
  read(40,*)
  if ((adjustl(trim(inspecies)).ne.'').and.(adjustl(trim(inspecies)).ne. &
       adjustl(trim(spsymb)))) goto 10
  write(*,'("Info(species): method: ",a,"; running Z = ",I4,", (",A,")")') &
 &  trim(apwdescr),nz, trim(spname)
  ! nuclear charge in units of e
  spzn=-dble(nz)
  ! minimum radial mesh point proportional to 1/sqrt(Z)
  sprmin=1.d-5
  ! set the number of radial mesh points proportional to number of nodes
  nrmt=200+50*(spn(spnst)-1)
  
  ! find the optimal effective infinity
  sprmax=80.d0
  radialgridtype="cubic"
   
  do i = 1, 2

! generate the radial mesh
    if ((radialgridtype.ne."cubic").and. &
        (radialgridtype.ne."expocubic").and. &
        (radialgridtype.ne."exponential")) then 
        write(*,*) 'Wrong radialGridType.'
        write(*,*) 'Choose between cubic, expocubic and exponential!'
        write(*,*) 'Terminating...'
        stop
    endif

! estimate the number of radial mesh points to infinity
    if (radialgridtype.eq."exponential") then
        t1 = Log (sprmax/sprmin) / Log (rmt/sprmin)
        t2 = dble (nrmt-1) * t1
        spnr = Nint (t2) + 1
    else
        spnr = 2+Nint(dble(nrmt-1)*((sprmax-sprmin)/(rmt-sprmin))**0.333333333333333d0)
    endif

! generate the radial meshes
    t1 = 1.d0 / dble (nrmt-1)
! logarithmic mesh
    t2 = Log (rmt/sprmin)
    cutoff=Nint(nrmt*0.5d0)

    if (allocated(r)) deallocate(r)
    allocate(r(spnr))
    Do ir = 1, spnr
        if (radialgridtype.eq."cubic") then
            r(ir) = sprmin+(dble(ir-1)/dble(nrmt-1))**3*(rmt-sprmin)
        elseif (radialgridtype.eq."exponential") then
            r(ir) = sprmin*Exp(dble(ir-1)*t1*t2)
        else
            r(ir) = 0.5d0*(erf(5d0*dble(ir-cutoff)/nrmt)+1d0)* &
           &  (sprmin+(dble(ir-1)/dble(nrmt-1))**3*(rmt-sprmin))+ &
           &  (1d0-0.5d0*(erf(5d0*dble(ir-cutoff)/nrmt)+1d0))*sprmin*Exp(dble(ir-1)*t1*t2)
        endif
    End Do

    if (allocated(rho)) deallocate(rho)
    if (allocated(vr)) deallocate(vr)
    if (allocated(rwf)) deallocate(rwf)
    allocate(rho(spnr))
    allocate(vr(spnr))
    allocate(rwf(spnr,2,spnst))

    ! solve the Kohn-Sham-Dirac equations for the atom
    call atom(.true., spzn, spnst, spn, spl, spk, spocc, xctype, xcgrad, &
    &         spnr, r, eval, rho, vr, rwf, nrmt, .true.)

    do ir=spnr,1,-1
      if (rho(ir).gt.1.d-10) then
        sprmax=1.5d0*r(ir)
        goto 20
      end if
    end do

  20 continue
  end do
  deallocate(rwf)
  
  ! check total charge is correct
  if (allocated(fr)) deallocate(fr)
  allocate(fr(spnr))
  if (allocated(gr)) deallocate(gr)
  allocate(gr(spnr))
  if (allocated(cf)) deallocate(cf)
  allocate(cf(3,spnr))
  do ir=1,spnr
    fr(ir)=4.d0*pi*rho(ir)*r(ir)**2
  end do
  call fderiv(-1,spnr,r,fr,gr,cf)
  if (abs(gr(spnr)+spzn).gt.1.d-5) then
    write(*,*)
    write(*,'("Error(species): charge mismatch")')
    write(*,*)
    stop
  end if
  deallocate(fr,gr,cf)
  
  ! find which states belong to core
  do ist=1,spnst
    if (eval(ist).lt.ecvcut) then
      spcore(ist)=.true.
    else
      spcore(ist)=.false.
    end if
  end do

  ! check that the state for same n and l but different k is also core
  do ist=1,spnst
    if (spcore(ist)) then
      do jst=1,spnst
        if ((spn(ist).eq.spn(jst)).and.(spl(ist).eq.spl(jst))) spcore(jst)=.true.
      end do
    end if
  end do

! default values of linearization energy in valence region
  elval(:) = 0.15d0

! Semi-core states  
  nl(:)=0
  maxl=0
  
  allocate(p0(nrmt),p1(nrmt),q0(nrmt),q1(nrmt))
  do ist = 1, spnst
    if (.not.spcore(ist)) then
      if ((spl(ist).eq.0).or.(spl(ist).eq.spk(ist))) then
        l=spl(ist)
        if (eval(ist).lt.esccut) then
        
          ! check whether the state is really SC one
          e = eval(ist)
          call rschroddme(0, l, 0, e, nrmt, r, vr, nn, p0, p1, q0, q1)
          
          t0 = p0(nrmt)
          
          do while (e<=elval(l))
            e = e+0.0001
            call rschroddme(0, l, 0, e, nrmt, r, vr, nn, p0, p1, q0, q1)
            
            t = p0(nrmt)
            if (t*t0<0.d0) then
              nl(l) = nl(l)+1
              el(l,nl(l)) = eval(ist)
              write(*,*) 'SEMI-CORE LOCAL ORBITAL: l=',l,'    e=', eval(ist)
              exit
            end if
            
          end do
        
        end if
      end if
      if (spl(ist).gt.maxl) maxl=spl(ist)
    end if
  end do
  maxl = min(maxl,3)
  deallocate(p0,p1,q0,q1)

!---------------------------------------------
! Advanced feature: LOs for unoccupied states
!---------------------------------------------

  if (locorbxs) then
    
    call findxslorb
    
    do l = 0, maxl
      do i = 1, nl(l)
        if (dabs(el(l,i))<=0.5d0) elval(l) = min( elval(l), el(l,i) )
      end do
    end do
    
  end if
  deallocate(vr)

!------------------------------
! Generate species.xml files
!------------------------------
  call writexmlspecies()

  ! read another element from file
  goto 10

end program
