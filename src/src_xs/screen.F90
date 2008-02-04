
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
  integer :: iq,ik
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
  call init0
  call init1
  call init2xs
  call readfermi
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)

  ! save variables for the Gamma q-point
  call tdsave0
  ! generate Gaunt coefficients
  call tdgauntgen(lmaxapw,lmaxemat,lmaxapw)
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)
  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)

  ! allocate arrays
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3,nst1,nst2))
  if (allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3,nst3,nst4))

  ! q-point parallelization
  call genparidxran('q')
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
        call ematqk1(iq,ik)
        call getpemat(iq,ik,trim(fnpmat),trim(fnemat),m12=xiou,m34=xiuo, &
             p12=pmou,p34=pmuo)

write(*,*) 'iq:',sum(abs(xiou)),sum(abs(xiuo)),sum(abs(pmou)),sum(abs(pmuo))

        ! accumulate matrix elements for dielectric matrix
        


        ! end loop over k-points
     end do



     ! store dielectric matrix for q-point


     ! end loop over q-points
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
