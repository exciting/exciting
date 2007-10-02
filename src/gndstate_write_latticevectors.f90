subroutine gndstate_write_latticevectors()
use modmain
use modmpi,only:rank
if(rank.eq.0) then
     call writelat
     ! write inter-atomic distances to file
     call writeiad
     ! write symmetry matrices to file
     call writesym
     ! output the k-point set to file
     call writekpts
     ! write lattice vectors and atomic positions to file
     call writegeom(.false.)
     ! open INFO.OUT file
     open(60,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
     ! open TOTENERGY.OUT
     open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
     ! open FERMIDOS.OUT
     open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
     ! open MOMENT.OUT if required
     if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
          form='FORMATTED')
     ! open FORCEMAX.OUT if required
     if (tforce) open(64,file='FORCEMAX'//trim(filext),action='WRITE', &
          form='FORMATTED')
     ! write out general information to INFO.OUT
     call writeinfo(60)
     ! initialise or read the charge density and potentials from file
     write(60,*)
	 call flushifc(60)
  endif
end subroutine