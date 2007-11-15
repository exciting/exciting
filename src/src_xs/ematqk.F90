
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ematqk
  implicit none
contains

  subroutine ematqk(iq,ik)
    use modmain
    use modxs
    use modmpi
    use m_tdintrg
    use m_ematqkgmt
    use m_ematqkgir
    use m_putemat
    use m_putdevalsv
    use m_emattim
    use m_getunit
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    ! local variables
    character(*), parameter :: thisnam = 'ematqk'
    character(256) :: filextt
    ! allocatable arrays
    complex(8), allocatable :: evecfvo(:,:),evecfvo0(:,:)
    complex(8), allocatable :: evecfvu(:,:),evecfvu0(:,:)
    complex(8), allocatable :: evecfvo2(:,:),evecfvo20(:,:)
    complex(8), allocatable :: evecfvu2(:,:),evecfvu20(:,:)
    complex(8), allocatable :: helpm(:,:),helpm2(:,:)
    integer :: ikq,igq, istc, istv
    integer :: i1,i2,i,j,recl,n,n0
    real(8) :: cpuini,cpuread,cpumain,cpuwrite,cpuall
    real(8) :: cpugnt,cpumt,cpuir
    real(8) :: cpumalores,cpumaloares,cpumloares,cpumloaares
    real(8) :: cpumlolores,cpumloloares,cpumirres,cpumirares,cpudbg
    real(8) :: cpu0,cpu1,cpu00,cpu01
    real(8) :: vql_(3), vkl_(3)
    complex(8) :: dum

    call cpu_time(cpu0)

    ikq=ikmapikq(iq,ik)

    ! check for stop statement
    write(msg,*) 'for q-point', iq, ': k-point:', ik-1, ' finished'
    call tdchkstop

    cpumtaa=0.d0; cpumtalo=0.d0; cpumtloa=0.d0; cpumtlolo=0.d0
    cpugnt=0.d0; cpumt=0.d0; cpuir=0.d0
    cpumalores=0.d0; cpumaloares=0.d0; cpumloares=0.d0; cpumloaares=0.d0
    cpumlolores=0.d0; cpumloloares=0.d0; cpumirres=0.d0; cpumirares=0.d0
    cpudbg=0.d0

    ! allocate temporary arrays
    n=ngk(ikq,1)
    n0=ngk0(ik,1)
    allocate(evecfvo0(nlotot,nstval))
    allocate(evecfvu0(nlotot,nstcon))
    allocate(evecfvo(nlotot,nstval))
    allocate(evecfvu(nlotot,nstcon))
    !
    allocate(evecfvo20(n0,nstval))
    allocate(evecfvu20(n0,nstcon))
    allocate(evecfvo2(n,nstval))
    allocate(evecfvu2(n,nstcon))
    !
    allocate(xihir(n0,n))
    allocate(helpm(nlotot,max(nstval,nstcon)))
    allocate(helpm2(n0,max(nstval,nstcon))) ! for ir

    ! read eigenvectors, eigenvalues and occupancies for G+k+q
    call getevecfv(vkl(1,ikq),vgkl(1,1,ikq,1),evecfv)
    call getevalsv(vkl(1,ikq),evalsv(1,ikq))
    ! assume non-metals here - occupation numbers equal for k and k+q
    call getoccsv(vkl(1,ikq),occsv(1,ikq))
    evecfvo(:,:) = evecfv(ngk(ikq,1)+1:ngk(ikq,1)+nlotot,1:nstval,1)
    evecfvu(:,:) = evecfv(ngk(ikq,1)+1:ngk(ikq,1)+nlotot,nstval+1:nstsv,1)
    evecfvo2(:,:) = evecfv(1:ngk(ikq,1),1:nstval,1)
    evecfvu2(:,:) = evecfv(1:ngk(ikq,1),nstval+1:nstsv,1)

    ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
    call genfilname(iq=0,setfilext=.true.)
    call getevecfv0(vkl0(1,ik),vgkl0(1,1,ik,1),evecfv0)
    call getevalsv0(vkl0(1,ik),evalsv0(1,ik))
    ! change back file extension
    call genfilname(iq=iq,setfilext=.true.)
    evecfvo0(:,:) = evecfv0(ngk0(ik,1)+1:ngk0(ik,1)+nlotot,1:nstval,1)
    evecfvu0(:,:) = evecfv0(ngk0(ik,1)+1:ngk0(ik,1)+nlotot,nstval+1:nstsv,1)
    evecfvo20(:,:) = evecfv0(1:ngk0(ik,1),1:nstval,1)
    evecfvu20(:,:) = evecfv0(1:ngk0(ik,1),nstval+1:nstsv,1)

    ! eigenvalue differences
    do istv=1,nstval
       do istc=1,nstcon
          ! resonant part
          deou(istv,istc) = evalsv0(istv,ik) - evalsv(nstval+istc,ikq)
          ! antiresonant part
          deuo(istc,istv) = evalsv0(nstval+istc,ik) - evalsv(istv,ikq)
       end do
    end do

    call cpu_time(cpu1)
    cpuini=cpu1-cpu0

    ! get expansion coefficients (q=0)
    call genfilname(basename='APWDLM',iq=0,filnam=fnevapw)
    inquire(iolength=recl) vql_,vkl_,apwdlm0
    call getunit(unit1)
    open(unit1,file=trim(fnevapw),action='read',&
         form='unformatted',status='old',access='direct',recl=recl)
    read(unit1,rec=ik) vql_,vkl_,apwdlm0
    close(unit1)
    ! get expansion coefficients (q)
    call genfilname(basename='APWDLM',iq=iq,filnam=fnevapw)
    inquire(iolength=recl) vql_,vkl_,apwdlm
    call getunit(unit1)
    open(unit1,file=trim(fnevapw),action='read',&
         form='unformatted',status='old',access='direct',recl=recl)
    read(unit1,rec=ikq) vql_,vkl_,apwdlm
    close(unit1)

    call cpu_time(cpu0)
    cpuread=cpu0-cpu1

    ! loop over G+q vectors
    do igq=1,ngq(iq)
       
       call terminate_inqr('ematqk')

       call cpu_time(cpu00)
       ! summation of Gaunt coefficients wrt radial integrals
       call tdintrg(iq,igq)
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
          !-------------------!
          !   resonant part   !
          !-------------------!
          ! multiplication xi = xiho * evecfvu
          call zgemm('n','n', nstval, nstcon, nlotot, zone, xiohalo, &
               nstval, evecfvu, nlotot, zone, xiou(1,1,igq), nstval )
          call cpu_time(cpu00)
          cpumalores=cpumalores+cpu00-cpu01
          !------------------------!
          !   anti-resonant part   !
          !------------------------!
          ! multiplication xi = xihu * evecfvo
          call zgemm('n','n', nstcon, nstval, nlotot, zone, xiuhalo, &
               nstcon, evecfvo, nlotot, zone, xiuo(1,1,igq), nstcon )
          call cpu_time(cpu01)
          cpumaloares=cpumaloares+cpu01-cpu00

          ! lo-APW contribution
          !-------------------!
          !   resonant part   !
          !-------------------!
          ! multiplication xi = evecfvo * xihu
          call zgemm('c','n', nstval, nstcon, nlotot, zone, evecfvo0, &
               nlotot, xiuhloa, nlotot, zone, xiou(1,1,igq), nstval )
          call cpu_time(cpu00)
          cpumloares=cpumloares+cpu00-cpu01
          !------------------------!
          !   anti-resonant part   !
          !------------------------!
          ! multiplication xi = evecfvu * xiho
          call zgemm('c','n', nstcon, nstval, nlotot, zone, evecfvu0, &
               nlotot, xiohloa, nlotot, zone, xiuo(1,1,igq), nstcon )
          call cpu_time(cpu01)
          cpumloaares=cpumloaares+cpu01-cpu00

          ! lo-lo contribution
          !-------------------!
          !   resonant part   !
          !-------------------!
          ! multiplication helpm(i,m) = xih * evecfvu
          call zgemm('n','n', nlotot, nstcon, nlotot, zone, xih, &
               nlotot, evecfvu, nlotot, zzero, helpm, nlotot )
          ! multiplication xi = hermc(evecfvo) * helpm
          call zgemm('c','n', nstval, nstcon, nlotot, zone, evecfvo0, &
               nlotot, helpm, nlotot, zone, xiou(1,1,igq), nstval )
          call cpu_time(cpu00)
          cpumlolores=cpumlolores+cpu00-cpu01
          !------------------------!
          !   anti-resonant part   !
          !------------------------!
          ! multiplication helpm(i,m) = xih * evecfvu
          call zgemm('n','n', nlotot, nstval, nlotot, zone, xih, &
               nlotot, evecfvo, nlotot, zzero, helpm, nlotot )
          ! multiplication xi = hermc(evecfvo) * helpm
          call zgemm('c','n', nstcon, nstval, nlotot, zone, evecfvu0, &
               nlotot, helpm, nlotot, zone, xiuo(1,1,igq), nstcon )
          call cpu_time(cpu01)
          cpumloloares=cpumloloares+cpu01-cpu00
       end if

       ! interstitial contribution
       !-------------------!
       !   resonant part   !
       !-------------------!
       ! multiplication helpm2(i,m) = xih * evecfvu2
       call zgemm('n','n', n0, nstcon, n, zone, xihir, &
            n0, evecfvu2, n, zzero, helpm2, n0 )
       ! multiplication xi = hermc(evecfvo2) * helpm2
       call zgemm('c','n', nstval, nstcon, n0, zone, evecfvo20, &
            n0, helpm2, n0, zone, xiou(1,1,igq), nstval )
       call cpu_time(cpu00)
       cpumirres=cpumirres+cpu00-cpu01
       !------------------------!
       !   anti-resonant part   !
       !------------------------!
       ! multiplication helpm2(i,m) = xih * evecfvu2
       call zgemm('n','n', n0, nstval, n, zone, xihir, &
            n0, evecfvo2, n, zzero, helpm2, n0 )
       ! multiplication xi = hermc(evecfvo2) * helpm2
       call zgemm('c','n', nstcon, nstval, n0, zone, evecfvu20, &
            n0, helpm2, n0, zone, xiuo(1,1,igq), nstcon )
       call cpu_time(cpu01)
       cpumirares=cpumirares+cpu01-cpu00

       if (dbglev.gt.0) then
          ! check non-diagonal parts of <phi_nk|exp(-i(G+q)r|phi_nk+q>
          do i1=1,nstval
             do i2=1,nstcon
                write(1000 + iq,'(4i6,3g18.10)') ik,igq,i1,i2, &
                     xiou(i1,i2,igq), &
                     abs(xiou(i1,i2,igq))**2
                write(2000 + iq,'(4i6,3g18.10)') ik,igq,i1,i2, &
                     xiuo(i2,i1,igq), &
                     abs(xiuo(i2,i1,igq))**2
             end do
          end do
       end if
       call cpu_time(cpu00)
       cpudbg=cpudbg+cpu00-cpu01
    end do ! igq
    call cpu_time(cpu1)
    cpumain=cpu1-cpu0


!************************************************************
! set anti-resonant term to zero for testing -> as switch in INPUT //////

!xiuo(:,:,:)=zzero

! simulate Lindhard function (set resonant term to one)
!xiou(:,:,:)=1.d0  

!************************************************************

    ! write to emat file
    call putemat(iq,ik,.false.,trim(fnemat_t),xiou,xiuo)

    ! write Kohn Sham energy differences
    call putdevalsv(iq,ik,.false.,trim(fndevalsv_t),deou,deuo)

    ! deallocate
    deallocate(helpm,xihir,evecfvo,evecfvu,evecfvo0,evecfvu0)
    deallocate(evecfvo2,evecfvu2,evecfvo20,evecfvu20)
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
  
end module m_ematqk
