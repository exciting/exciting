
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
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
  use m_getpemat
  implicit none
  ! local variables
  character(*), parameter :: thisnam='screen'
  real(8), allocatable :: scis12(:,:),scis21(:,:)
  real(8) :: vklofft(3),rgkmaxt
  integer :: iq,ik,ikq,ngridkt(3),nemptyt,nwdft
  logical :: nosymt,reducekt
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  rgkmaxt=rgkmax
  nemptyt=nempty
  nwdft=nwdf
  ! map variables for screening
  call initscr
  nosym=nosymscr
  ! no symmetries implemented for screening
  reducek=.false.
  ngridk(:)=ngridkscr(:)
  vkloff(:)=vkloffscr(:)
  rgkmax=rgkmaxscr
  nempty=nemptyscr
  ! only one frequency w=0
  nwdf=1
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)

  ! call dielectric function with only one frequency point
  call df

!!$  ! typ of matrix elements
!!$  emattype=1
!!$  call init0
!!$  call init1
!!$  call init2xs
!!$  call readfermi
!!$  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
!!$  call genfilname(basename='PMAT',appfilext=.true.,filnam=fnpmat)
!!$  ! save variables for the Gamma q-point
!!$  call tdsave0
!!$  ! generate Gaunt coefficients
!!$  call tdgauntgen(lmaxapw,lmaxemat,lmaxapw)
!!$  ! find indices for non-zero Gaunt coefficients
!!$  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)
!!$  ! find highest (partially) occupied and lowest (partially) unoccupied states
!!$  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
!!$  call ematbdcmbs(emattype)
!!$  ! allocate arrays
!!$  if (allocated(deou)) deallocate(deou)
!!$  allocate(deou(nst1,nst2))
!!$  if(allocated(deuo)) deallocate(deuo)
!!$  allocate(deuo(nst3,nst4))
!!$  if (allocated(docc12)) deallocate(docc12)
!!$  allocate(docc12(nst1,nst2))
!!$  if (allocated(docc21)) deallocate(docc21)
!!$  allocate(docc21(nst3,nst4))
!!$  if (allocated(pmou)) deallocate(pmou)
!!$  allocate(pmou(3,nst1,nst2))
!!$  if (allocated(pmuo)) deallocate(pmuo)
!!$  allocate(pmuo(3,nst3,nst4))
!!$  allocate(scis12(nst1,nst2))
!!$  allocate(scis21(nst2,nst1))
!!$  ! q-point parallelization
!!$  call genparidxran('q')
!!$  do iq=qpari,qparf
!!$     write(*,*) 'q-loop:iq',iq
!!$     call genfilname(nodotpar=.true.,basename='EMAT',iq=iq,&
!!$          etype=emattype,procs=procs,rank=rank,appfilext=.true.,filnam=fnetim)
!!$     call updateq(iq)
!!$     ! calculate k+q and G+k+q related variables
!!$     call init1xs(qvkloff(1,iq))
!!$     ! write G+q-vectors
!!$     call writegqpts(iq)
!!$     call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
!!$     ! generate radial integrals wrt. sph. Bessel functions
!!$     call ematrad(iq)
!!$     call ematqalloc
!!$     do ik=1,nkpt
!!$        ikq=ikmapikq(ik,iq)
!!$        call getdevaldoccsv(iq,ik,ikq,istlo1,isthi1,istlo2,isthi2,deou,docc12, &
!!$             scis12)
!!$        call getdevaldoccsv(iq,ik,ikq,istlo2,isthi2,istlo1,isthi1,deuo,docc21, &
!!$             scis21)
!!$        call ematqk1(iq,ik)
!!$        call getpemat(iq,ik,trim(fnpmat),trim(fnemat),m12=xiou,m34=xiuo, &
!!$             p12=pmou,p34=pmuo)
!!$
!!$write(*,'(a,i6,4f12.4)') 'ik:',ik,sum(abs(xiou)),sum(abs(xiuo)),sum(abs(pmou)),sum(abs(pmuo))
!!$
!!$        ! accumulate matrix elements for dielectric matrix
!!$        ! end loop over k-points
!!$     end do
!!$     call ematqdealloc
!!$     ! store dielectric matrix for q-point
!!$     ! end loop over q-points
!!$  end do
!!$  deallocate(xiou,xiuo,pmou,pmuo,deou,deuo,docc12,docc21,scis12,scis21)

  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  rgkmax=rgkmaxt
  nempty=nemptyt
  nwdf=nwdft
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine screen
