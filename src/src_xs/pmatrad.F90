 
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pmatrad
  use modmain
  use modxs
  use m_gradzfmtr
  use m_getunit
  implicit none
  ! local variables
  integer :: is,ia,ias,nr,ir
  integer :: l1,m1,lm1,l3,m3,lm3
  integer :: ilo,ilo1,ilo2,io,io1,io2,j
  real(8) :: cpu0,cpu1
  integer :: u11,u22,u33,u44
  ! automatic arrays
  real(8) :: r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
  ! allocatable arrays
  real(8), allocatable :: fapw(:),flo(:),dapwfr(:,:,:,:,:),dlofr(:,:,:,:,:)

  ! allocate local arrays for radial derivatives
  allocate(fapw(nrmtmax))
  allocate(flo(nrmtmax))
  allocate(dapwfr(lmmaxapw,nrmtmax,3,apwordmax,lmmaxapw))
  allocate(dlofr(lmmaxapw,nrmtmax,3,nlomax,-lolmax:lolmax))
  ! zero arrays for radial derivatives
  dapwfr(:,:,:,:,:)=0.d0
  dlofr(:,:,:,:,:)=0.d0

  ! allocate arrays for radial integrals
  if (allocated(ripaa)) deallocate(ripaa)
  if (allocated(ripalo)) deallocate(ripalo)
  if (allocated(riploa)) deallocate(riploa)
  if (allocated(riplolo)) deallocate(riplolo)
  allocate(ripaa(apwordmax,lmmaxapw,apwordmax,lmmaxapw,natmtot,3))
  allocate(ripalo(apwordmax,lmmaxapw,nlomax,-lolmax:lolmax,natmtot,3))
  allocate(riploa(nlomax,-lolmax:lolmax,apwordmax,lmmaxapw,natmtot,3))
  allocate(riplolo(nlomax,-lolmax:lolmax,nlomax,-lolmax:lolmax,natmtot,3))
  ! zero arrays for radial integrals
  ripaa(:,:,:,:,:,:)=0.d0
  ripalo(:,:,:,:,:,:)=0.d0
  riploa(:,:,:,:,:,:)=0.d0
  riplolo(:,:,:,:,:,:)=0.d0

  if (dbglev.gt.1) then
     ! APW-APW
     call getunit(u11)
     open(unit=u11,file='IRADPaa'//filext,form='formatted',action='write', &
          status='replace')
     write(u11,'(a)') 'igq,ias,l1,io1,l3,io2,l2   iraa'
     write(u11,'(a)') '-----------------------------------------------------'
     ! APW-lo
     call getunit(u22)
     open(unit=u22,file='IRADPalo'//filext,form='formatted',action='write', &
          status='replace')
     write(u22,'(a)') 'igq,ias,ilo,l1,l3,io,l2,   iralo'
     write(u22,'(a)') '-----------------------------------------------------'
     ! lo-APW
     call getunit(u33)
     open(unit=u33,file='IRADPloa'//filext,form='formatted',action='write', &
          status='replace')
     write(u33,'(a)') 'igq,ias,ilo,l1,l3,io,l2,   irloa'
     write(u33,'(a)') '-----------------------------------------------------'
     ! lo-lo
     call getunit(u44)
     open(unit=u44,file='IRADPlolo'//filext,form='formatted',action='write', &
          status='replace')
     write(u44,'(a)') 'igq,ias,ilo1,l1,ilo2,l3,l2,   irlolo'
     write(u44,'(a)') '-----------------------------------------------------'
  end if

  ! begin loop over species
  do is=1,nspecies
     nr=nrmt(is)
     do ir=1,nr
        ! calculate r^2
        r2(ir)=spr(ir,is)**2
     end do
     ! begin loop over atoms
     do ia=1,natoms(is)
        ias=idxas(ia,is)

        call cpu_time(cpu0)

        !--------------------!
        !     derivatives    !
        !--------------------!
        ! APW functions
        do l1=0,lmaxapw
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              do io=1,apword(l1,is)
                 fapw(:)=apwfr(:,1,io,l1,ias)
                 call gradzfmtr(lmaxapw,nr,spr(1,is),l1,m1,lmmaxapw, &
                      nrmtmax,fapw,dapwfr(1,1,1,io,lm1))
             end do
           end do
        end do
        ! local orbital functions
        do ilo=1,nlorb(is)
           l1=lorbl(ilo,is)
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              flo(:)=lofr(:,1,ilo,ias)
              call gradzfmtr(lmaxapw,nr,spr(1,is),l1,m1,lmmaxapw, &
                   nrmtmax,flo,dlofr(1,1,1,ilo,m1))
           end do
        end do

call cpu_time(cpu1)
write(*,'(a,i6,f12.3)') 'deriv: ',ias,cpu1-cpu0
!!$write(555,*) dapwfr
!!$write(556,*) dlofr

        !----------------!
        !     APW-APW    !
        !----------------!
        do l1=0,lmaxapw
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              do io1=1,apword(l1,is)
                 do l3=0,lmaxapw
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       do io2=1,apword(l3,is)
                          do j=1,3
                             fr(:)=r2(1:nr)*apwfr(1:nr,1,io1,l1,ias)* &
                                  dapwfr(lm1,1:nr,j,io2,lm3)
                             call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                             ripaa(io1,lm1,io2,lm3,ias,j)=gr(nr)
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do

call cpu_time(cpu0)
write(*,'(a,i6,f12.3)') 'AA  : ',ias,cpu0-cpu1
if (ias.eq.1) then
do lm1=1,2
do lm3=1,2
   write(557,'(2i6,3g18.10)') lm1,lm3,ripaa(1,lm1,1,lm3,1,:)
end do
end do
end if

        !----------------------------!
        !     APW-local-orbital      !
        !----------------------------!
        do l1=0,lmaxapw
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              do io1=1,apword(l1,is)
                 do ilo=1,nlorb(is)
                    l3=lorbl(ilo,is)
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       do j=1,3
                          fr(:)=r2(1:nr)*apwfr(:,1,io1,l1,ias)* &
                               dlofr(lm1,1:nr,j,ilo,m3)
                          call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                          ripalo(io1,lm1,ilo,m3,ias,j)=gr(nr)
                       end do
                    end do
                 end do
              end do
           end do
        end do

call cpu_time(cpu1)
write(*,'(a,i6,f12.3)') 'Alo : ',ias,cpu1-cpu0
write(558,*) ripalo

        !----------------------------!
        !     local-orbital-APW      !
        !----------------------------!
        do ilo=1,nlorb(is)
           l1=lorbl(ilo,is)
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              do l3=0,lmaxapw
                 do m3=-l3,l3
                    lm3=idxlm(l3,m3)
                    do io2=1,apword(l3,is)
                       do j=1,3
                          fr(:)=r2(1:nr)*lofr(:,1,ilo,ias)* &
                               dapwfr(lm1,1:nr,j,io2,lm3)
                          call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                          riploa(ilo,m1,io2,lm3,ias,j)=gr(nr)
                       end do
                    end do
                 end do
              end do
           end do
        end do

call cpu_time(cpu0)
write(*,'(a,i6,f12.3)') 'loA : ',ias,cpu0-cpu1
write(559,*) riploa

        !------------------------------------!
        !     local-orbital-local-orbital    !
        !------------------------------------!
        do ilo1=1,nlorb(is)
           l1=lorbl(ilo1,is)
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              do ilo2=1,nlorb(is)
                 l3=lorbl(ilo2,is)
                 do m3=-l3,l3
                    lm3=idxlm(l3,m3)
                    do j=1,3
                       fr(:)=r2(1:nr)*lofr(:,1,ilo1,ias)* &
                            dlofr(lm1,1:nr,j,ilo2,m3)
                       call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                       riplolo(ilo1,m1,ilo2,m3,ias,j)=gr(nr)
                    end do
                 end do
              end do
           end do
        end do

call cpu_time(cpu1)
write(*,'(a,i6,f12.3)') 'lolo: ',ias,cpu1-cpu0
write(560,*) riplolo

        ! end loops over atoms and species
     end do
  end do

  ! deallocate
  deallocate(fapw,flo,dapwfr,dlofr)
  if (dbglev.gt.1) then
     ! close files
     close(u11)
     close(u22)
     close(u33)
     close(u44)
  end if

end subroutine pmatrad
