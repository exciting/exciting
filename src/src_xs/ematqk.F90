
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqk(iq,ik)
  use modmain
  use modxs
  use modmpi
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
  complex(8), allocatable :: helpm(:,:),helpm2(:,:)
  integer :: ikq,igq
  integer :: i1,i2,recl,n,n0
  real(8) :: cpuini,cpuread,cpumain,cpuwrite,cpuall
  real(8) :: cpugnt,cpumt,cpuir
  real(8) :: cpumalores,cpumaloares,cpumloares,cpumloaares
  real(8) :: cpumlolores,cpumloloares,cpumirres,cpumirares,cpudbg
  real(8) :: cpu0,cpu1,cpu00,cpu01
  real(8) :: vql_(3), vkl_(3)

  call cpu_time(cpu0)
  ! find equivalent k-point
  ikq=ikmapikq(ik,iq)
  !    if ((modulo(ik,nkpt/10).eq.0).or.(ik.eq.nkpt)) &
  !         write(*,'("Info(ematqk2): ",I6,I6," of ",I6," k-points")') ik,ikq,nkpt
  ! check for stop statement
  write(msg,*) 'for q-point', iq, ': k-point:', ik-1, ' finished'
  call tdchkstop

  cpumtaa=0.d0; cpumtalo=0.d0; cpumtloa=0.d0; cpumtlolo=0.d0
  cpugnt=0.d0; cpumt=0.d0; cpuir=0.d0
  cpumalores=0.d0; cpumaloares=0.d0; cpumloares=0.d0; cpumloaares=0.d0
  cpumlolores=0.d0; cpumloloares=0.d0; cpumirres=0.d0; cpumirares=0.d0
  cpudbg=0.d0

  ! allocate temporary arrays
  n0=ngk0(ik,1)
  n=ngk(ikq,1)
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
  allocate(helpm(nlotot,max(nst1,nst2)))
  allocate(helpm2(n0,max(nst1,nst2)))
  ! zero arrays
  xiohalo(:,:)=zzero
  xiuhloa(:,:)=zzero

  ! read eigenvectors, eigenvalues and occupancies for G+k+q
  call getevecfv(vkl(1,ikq),vgkl(1,1,ikq,1),evecfv)
  call getevalsv(vkl(1,ikq),evalsv(1,ikq))
  ! read occupation numbers for G+k+q
  call getoccsv(vkl(1,ikq),occsv(1,ikq))
  evecfvu(:,:)=evecfv(ngk(ikq,1)+1:ngk(ikq,1)+nlotot,istlo2:isthi2,1)
  evecfvu2(:,:)=evecfv(1:ngk(ikq,1),istlo2:isthi2,1)

  ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
!!!!!!!!!!!!  call genfilname(iqmt=0,setfilext=.true.)
  call getevecfv0(vkl0(1,ik),vgkl0(1,1,ik,1),evecfv0)
  call getevalsv0(vkl0(1,ik),evalsv0(1,ik))
  ! read occupation numbers for G+k
  call getoccsv0(vkl0(1,ik),occsv0(1,ik))
  evecfvo0(:,:)=evecfv0(ngk0(ik,1)+1:ngk0(ik,1)+nlotot,istlo1:isthi1,1)
  evecfvo20(:,:)=evecfv0(1:ngk0(ik,1),istlo1:isthi1,1)
  ! change back file extension
!!!!!!!!!!!  call genfilname(iqmt=iq,setfilext=.true.)

  call cpu_time(cpu1)
  cpuini=cpu1-cpu0

!!$  ! get expansion coefficients (q=0)
!!$  call genfilname(basename='APWDLM',iqmt=0,filnam=fnevapw)
!!$  inquire(iolength=recl) vql_,vkl_,apwdlm0
!!$  call getunit(unit1)
!!$  open(unit1,file=trim(fnevapw),action='read',&
!!$       form='unformatted',status='old',access='direct',recl=recl)
!!$  read(unit1,rec=ik) vql_,vkl_,apwdlm0
!!$  close(unit1)
!!$  ! get expansion coefficients (q)
!!$  call genfilname(basename='APWDLM',iqmt=iq,filnam=fnevapw)
!!$  inquire(iolength=recl) vql_,vkl_,apwdlm
!!$  call getunit(unit1)
!!$  open(unit1,file=trim(fnevapw),action='read',&
!!$       form='unformatted',status='old',access='direct',recl=recl)
!!$  read(unit1,rec=ikq) vql_,vkl_,apwdlm
!!$  close(unit1)

  call cpu_time(cpu0)
  cpuread=cpu0-cpu1

  ! zero matrix elements array
  xiou(:,:,:)=zzero

  ! loop over G+q vectors
  do igq=1,ngq(iq)
     call terminate_inqr('ematqk')

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

     if (nlotot.gt.0) then 
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
        call cpu_time(cpu00)
        cpumloares=cpumloares+cpu00-cpu01

        ! lo-lo contribution
        ! multiplication helpm(i,m) = xih * evecfvu
        call zgemm('n','n', nlotot, nst2, nlotot, zone, xih, &
             nlotot, evecfvu, nlotot, zzero, helpm, nlotot )
        ! multiplication xi = hermc(evecfvo) * helpm
        call zgemm('c','n', nst1, nst2, nlotot, zone, evecfvo0, &
             nlotot, helpm, nlotot, zone, xiou(1,1,igq), nst1 )
        call cpu_time(cpu00)
        cpumlolores=cpumlolores+cpu00-cpu01
     end if

     ! interstitial contribution
     ! multiplication helpm2(i,m) = xih * evecfvu2
     call zgemm('n','n', n0, nst2, n, zone, xihir, &
          n0, evecfvu2, n, zzero, helpm2, n0 )
     ! multiplication xi = hermc(evecfvo2) * helpm2
     call zgemm('c','n', nst1, nst2, n0, zone, evecfvo20, &
          n0, helpm2, n0, zone, xiou(1,1,igq), nst1 )
     call cpu_time(cpu00)
     cpumirres=cpumirres+cpu00-cpu01

     if (dbglev.gt.0) then
        ! check non-diagonal parts of <phi_nk|exp(-i(G+q)r|phi_nk+q>
        do i1=1,nst1
           do i2=1,nst2
              write(1100 + iq,'(a,4i6,3g18.10)') 'ik,igq,i1,i2', &
                   ik,igq,i1,i2,xiou(i1,i2,igq),abs(xiou(i1,i2,igq))**2
           end do
        end do
     end if
     call cpu_time(cpu00)
     cpudbg=cpudbg+cpu00-cpu01
  end do ! igq

  call cpu_time(cpu1)
  cpumain=cpu1-cpu0

  ! deallocate
  deallocate(helpm,helpm2,xihir)
  deallocate(evecfvu,evecfvo0)
  deallocate(evecfvu2,evecfvo20)
  call cpu_time(cpu0)
  cpuwrite=cpu0-cpu1
  cpuall=cpuini+cpuread+cpumain+cpuwrite

  ! write timing information
  call emattim(iq,ik,trim(fnetim),&
       cpuini,cpuread,cpumain,cpuwrite,cpuall, &
       cpugnt,cpumt,cpuir, &
       cpumalores,cpumaloares,cpumloares,cpumloaares, &
       cpumlolores,cpumloloares,cpumirres,cpumirares,cpudbg, &
       cpumtaa,cpumtalo,cpumtloa,cpumtlolo)

end subroutine ematqk
