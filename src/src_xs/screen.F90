
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine screen
  use modmain
  use modmpi
  use modxs
  use m_filedel
  use m_genfilname
  use m_tdgauntgen
  use m_findgntn0
  use m_writegqpts
  implicit none
  ! local variables
  character(*), parameter :: thisnam='screen'
  integer :: taskt,iq,ik
  real(8) :: vklofft(3),rgkmaxt
  integer :: ngridkt(3),nemptyt
  logical :: nosymt,reducekt
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  rgkmaxt=rgkmax
  nemptyt=nempty
  ! map variables for screening
  call initscr
  nosym=nosymscr
  ! no symmetries implemented for screening
  reducek=.false.
  ngridk(:)=ngridkscr(:)
  vkloff(:)=vkloffscr(:)
  rgkmax=rgkmaxscr
  nempty=nemptyscr
  ! only one SCF iteration
  maxscl=1
  call init0
  call init1
  call init2xs
  taskt=task
  task=1
  isreadstate0=.true.
  call genfilname(setfilext=.true.,dotext='_SCR.OUT')
  ! calculate eigenvectors, -values and occupancies for basic k-mesh
  call gndstate
  task=taskt
  if (rank.eq.0) then
     ! safely remove unnecessary files
     call filedel('EQATOMS'//trim(filext))
     call filedel('EVALCORE'//trim(filext))
     call filedel('FERMIDOS'//trim(filext))
     call filedel('GEOMETRY'//trim(filext))
     call filedel('LATTICE'//trim(filext))
     call filedel('IADIST'//trim(filext))
     call filedel('LINENGY'//trim(filext))
     call filedel('SYMCRYS'//trim(filext))
     call filedel('SYMLAT'//trim(filext))
     call filedel('SYMSITE'//trim(filext))
     call filedel('TOTENERGY'//trim(filext))
     call filedel('EVALFV'//trim(filext))
     call filedel('RMSDVEFF'//trim(filext))
  end if

  partype='q'
  call genparidxran(partype)
  if (rank.eq.0) call writeqpts
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call tdsave0
  ! generate Gaunt coefficients
  call tdgauntgen(lmaxapw,lmaxemat,lmaxapw)
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)
  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)

  do iq=qpari,qparf
write(*,*) 'q-loop:iq',iq
     call updateq(iq)
     ! calculate k+q and G+k+q related variables
     call init1xs(qvkloff(1,iq))
     ! write G+q-vectors
     call writegqpts(iq)
     call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     ! generate radial integrals wrt. sph. Bessel functions
     call ematrad(iq)

     do ik=1,nkpt
     ! generate matrix elements of plane wave and momentum operator


     ! accumulate matrix elements for dielectric matrix
     end do



     ! store dielectric matrix for q-point


  end do

write(*,*) '------------------------------------------------------'

  ! *** TEST ***
  do iq=1,nqpt
     call updateq(iq)
     write(*,'(a,i6,3f12.3,3x,3f12.3)') 'TEST: iq/vql/vqlcu',iq,vql(:,iq),vqlcu
     call init1
  end do

  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  rgkmax=rgkmaxt
  nempty=nemptyt
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine screen
