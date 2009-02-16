
! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqk(iq,ik)
  use modmain
  use modmpi
  use modxs
  use summations
  use m_getapwcmt
  use m_getlocmt
  use m_putemat
  use m_emattim
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq,ik
  ! local variables
  character(*), parameter :: thisnam='ematqk'
  ! allocatable arrays
  complex(8), allocatable :: evecfvo0(:,:)
  complex(8), allocatable :: evecfvu(:,:)
  complex(8), allocatable :: evecfvo20(:,:)
  complex(8), allocatable :: evecfvu2(:,:)
  integer :: ikq,igq,n,n0
  real(8) :: cpuini,cpuread,cpumain,cpuwrite,cpuall
  real(8) :: cpugnt,cpumt,cpuir
  real(8) :: cpumalores,cpumloares
  real(8) :: cpumlolores,cpumirres,cpudbg
  real(8) :: cpu0,cpu1,cpu00,cpu01

  if (task.eq.330) call chkpt(3,(/task,iq,ik/),'ematqk: task, q-point index, &
       &k-point index; q-dependent matrix elements')
  call cpu_time(cpu0)
  ! find k+q-point
  ikq=ikmapikq(ik,iq)
  ! check for stop statement
  write(msg,*) 'for q-point', iq, ': k-point:', ik-1, ' finished'
  call xschkstop

  cpumtaa=0.d0; cpumtalo=0.d0; cpumtloa=0.d0; cpumtlolo=0.d0
  cpugnt=0.d0; cpumt=0.d0; cpuir=0.d0
  cpumalores=0.d0; cpumloares=0.d0
  cpumlolores=0.d0; cpumirres=0.d0
  cpudbg=0.d0

  ! allocate temporary arrays
  n0=ngk0(1,ik)
  n=ngk(1,ikq)
  ! allocate matrix elements array
  if (allocated(xiohalo)) deallocate(xiohalo)
  allocate(xiohalo(nst1,nlotot))
  if (allocated(xiuhloa)) deallocate(xiuhloa)
  allocate(xiuhloa(nlotot,nst2))
  ! allocate temporary arrays
  allocate(evecfvo0(nlotot,nst1))
  allocate(evecfvu(nlotot,nst2))
  allocate(evecfvo20(n0,nst1))
  allocate(evecfvu2(n,nst2))
  allocate(xihir(n0,n))
  ! zero arrays
  xiohalo(:,:)=zzero
  xiuhloa(:,:)=zzero

  call cpu_time(cpu1)
  cpuini=cpu1-cpu0

  ! read eigenvectors, eigenvalues and occupancies for G+k+q
  call getevecfv(vkl(1,ikq),vgkl(1,1,1,ikq),evecfv)
  call getevalsv(vkl(1,ikq),evalsv(1,ikq))
  ! read occupation numbers for G+k+q
  call getoccsv(vkl(1,ikq),occsv(1,ikq))
  evecfvu(:,:)=evecfv(ngk(1,ikq)+1:ngk(1,ikq)+nlotot,istl2:istu2,1)
  evecfvu2(:,:)=evecfv(1:ngk(1,ikq),istl2:istu2,1)

  ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
  call getevecfv0(vkl0(1,ik),vgkl0(1,1,1,ik),evecfv0)
  call getevalsv0(vkl0(1,ik),evalsv0(1,ik))
  ! read occupation numbers for G+k
  call getoccsv0(vkl0(1,ik),occsv0(1,ik))
  evecfvo0(:,:)=evecfv0(ngk0(1,ik)+1:ngk0(1,ik)+nlotot,istl1:istu1,1)
  evecfvo20(:,:)=evecfv0(1:ngk0(1,ik),istl1:istu1,1)
  ! change back file extension

  call getapwcmt(0,ik,1,nstfv,lmaxapwwf,apwcmt0)
  call getapwcmt(iq,ikq,1,nstfv,lmaxapwwf,apwcmt)
  call getlocmt(0,ik,1,nstfv,locmt0)
  call getlocmt(iq,ikq,1,nstfv,locmt)

  call cpu_time(cpu0)
  cpuread=cpu0-cpu1

  ! zero matrix elements array
  xiou(:,:,:)=zzero

  ! loop over G+q vectors
  do igq=1,ngq(iq)
     call terminateqry('ematqk')
     call cpu_time(cpu00)
     ! summation of Gaunt coefficients wrt radial integrals
     call ematgntsum(iq,igq)
     call cpu_time(cpu01)
     cpugnt=cpugnt+cpu01-cpu00
     ! muffin-tin contribution
     call ematqkgmt(iq,ik,igq)
     call cpu_time(cpu00)
     cpumt=cpumt+cpu00-cpu01
     ! interstitial contribution
     call ematqkgir(iq,ik,igq)
     call cpu_time(cpu01)
     cpuir=cpuir+cpu01-cpu00

     if ((.not.fastemat).and.(nlotot.gt.0)) then 
        ! muffin-tin contributions
        ! APW-lo contribution
        ! multiplication xi = xiho * evecfvu
        call zgemm('n','n', nst1, nst2, nlotot, zone, xiohalo, &
             nst1, evecfvu, nlotot, zone, xiou(1,1,igq), nst1 )
        call cpu_time(cpu00)
        cpumalores=cpumalores+cpu00-cpu01
        ! lo-APW contribution
        ! multiplication xi = evecfvo * xihu
        call zgemm('c','n', nst1, nst2, nlotot, zone, evecfvo0, &
             nlotot, xiuhloa, nlotot, zone, xiou(1,1,igq), nst1 )
        call cpu_time(cpu01)
        cpumloares=cpumloares+cpu01-cpu00
        ! lo-lo contribution
        call doublesummation_simple_cz(xiou(:,:,igq),evecfvo0,xih,evecfvu, &
             zone,zone,.true.)

        call cpu_time(cpu00)
        cpumlolores=cpumlolores+cpu00-cpu01
        cpu01=cpu00
     end if

     ! interstitial contribution
     call doublesummation_simple_cz(xiou(:,:,igq),evecfvo20,xihir,evecfvu2, &
          zone,zone,.true.)

     call cpu_time(cpu00)
     cpumirres=cpumirres+cpu00-cpu01

     call cpu_time(cpu01)
     cpudbg=cpudbg+cpu01-cpu00
  end do ! igq

  call cpu_time(cpu1)
  cpumain=cpu1-cpu0

  ! deallocate
  deallocate(xihir)
  deallocate(evecfvu,evecfvo0)
  deallocate(evecfvu2,evecfvo20)
  call cpu_time(cpu0)
  cpuwrite=cpu0-cpu1
  cpuall=cpuini+cpuread+cpumain+cpuwrite

  ! write timing information
  if ((task.ne.430).and.(task.ne.440).and.(task.ne.441).and.(task.ne.450) &
       .and.(task.ne.451)) then
     call emattim(iq,ik,trim(fnetim),cpuini,cpuread,cpumain,cpuwrite,cpuall, &
          cpugnt,cpumt,cpuir,cpumalores,cpumloares,cpumlolores,cpumirres, &
          cpudbg,cpumtaa,cpumtalo,cpumtloa,cpumtlolo)
  end if

end subroutine ematqk
