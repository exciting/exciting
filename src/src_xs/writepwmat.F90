
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepwmat
! !INTERFACE:
subroutine writepwmat
  ! !USES:
  use modmain
  use modxs
  ! !DESCRIPTION:
  !   Calculates the matrix elements of the plane wave $exp(-i(\mathbf{G}+
  !   mathbf{q})\mathbf{r})$ using routine {\tt genpwmat} and writes them to
  !   direct access file {\tt PWMAT.OUT}.
  !
  ! !REVISION HISTORY:
  !   Created November 2007 (Sagmeister)
  !EOP
  !BOC
  implicit none
  ! local variables
  integer, parameter :: iq=1
  integer ik,ikp,recl,isymkp,igq,un,un2
  real(8) :: vpl(3),vkpl(3)
  complex(8), allocatable :: apwalmk(:,:,:,:),apwalmkp(:,:,:,:)
  complex(8), allocatable :: evecfvk(:,:),evecfvkp(:,:)
  complex(8), allocatable :: evecsvk(:,:),evecsvkp(:,:)
  complex(8), allocatable :: pwmat(:,:,:),pwmatf(:,:,:)

    integer :: nstval_, nstcon_, nkpt_, ngq_,reclfull
    real(8) :: vql_(3), vkl_(3)

  integer :: j,iknr,lspl,isym,jsym,s(3,3),vg(3),ig,ist,jst, si(3,3)
  real(8) :: c(3,3),vt(3),v1(3),t1,t2,t3
  complex(8) :: zt1,zt2
  complex(8), allocatable :: yiou(:,:,:),yiuo(:,:,:)
  real(8), parameter :: epsrot=1.d-12

  ! initialise universal variables
  call init0
  call init1
  call init2xs
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  allocate(apwalmk(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(apwalmkp(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(evecfvk(nmatmax,nstfv))
  allocate(evecfvkp(nmatmax,nstfv))
  allocate(evecsvk(nstsv,nstsv))
  allocate(evecsvkp(nstsv,nstsv))
  ! allocate the momentum matrix elements array
  allocate(pwmat(ngq(iq),nstsv,nstsv))
  allocate(pwmatf(ngq(iq),nstsv,nstsv))

  ! allocate matrix elements array
  if (allocated(xiou)) deallocate(xiou)
  allocate(xiou(nstval,nstcon,ngq(iq)))
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiuo(nstcon,nstval,ngq(iq)))

  if (allocated(yiou)) deallocate(yiou)
  allocate(yiou(nstval,nstcon,ngq(iq)))
  if (allocated(yiuo)) deallocate(yiuo)
  allocate(yiuo(nstcon,nstval,ngq(iq)))

  ! read in the density and potentials from file
  call readstate
  ! find the new linearisation energies
  call linengy
  ! generate the APW radial functions
  call genapwfr
  ! generate the local-orbital radial functions
  call genlofr
  ! find the record length
  inquire(iolength=recl) pwmat
!!$  open(50,file='PWMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
!!$       status='REPLACE',recl=recl)
  open(50,file='PWMAT_ASC.OUT',action='WRITE',form='FORMATTED',status='REPLACE')
  do ik=1,nkpt
     write(*,*) 'Info(writepwmat): ik',ik
     ! get the eigenvectors from file for k-point
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvk)
     call getevecsv(vkl(1,ik),evecsvk)
     ! find the matching coefficients for k-point
     call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalmk)
     ! find k-point equivalent to k+q
     vpl(:)=vql(:,iq)+vkl(:,ik)
     call findkpt(vpl,isymkp,ikp)
     write(*,*) 'ik,ikp',ik,ikp
     vkpl(:)=vkl(:,ikp)
     ! get the eigenvectors from file for kp-point
     call getevecfv(vkl(1,ikp),vgkl(1,1,ikp,1),evecfvkp)
     call getevecsv(vkl(1,ikp),evecsvkp)
     ! find the matching coefficients for kp-point
     call match(ngk(ikp,1),gkc(1,ikp,1),tpgkc(1,1,ikp,1),sfacgk(1,1,ikp,1), &
          apwalmkp)
     ! calculate the matrix elements of the plane wave
     call genpwmat(vql(1,iq),ngqmax,ngq(iq),vgqc(1,1,iq),gqc(1,iq),igqig(1,iq),&
          ylmgq(1,1,iq),sfacgq(1,1,iq),vkl(1,ik),ngk(ik,1),igkig(1,ik,1), &
          apwalmk,evecfvk,evecsvk,vkl(1,ikp),ngk(ikp,1),igkig(1,ikp,1), &
          apwalmkp,evecfvkp,evecsvkp,pwmat)
     ! write to ASCII file
     do igq=1,ngq(iq)
        do ist=1,nstsv
           do jst=1,nstsv
              write(50,'(4i5,3g18.10)') ik,igq,ist,jst,pwmat(igq,ist,jst), &
                   abs(pwmat(igq,ist,jst))**2
           end do
        end do
     end do
     do igq=1,ngq(iq)
        xiou(:,:,igq)=pwmat(igq,1:nstval,nstval+1:nstsv)
        xiuo(:,:,igq)=pwmat(igq,nstval+1:nstsv,1:nstval)
     end do
     ! write to direct access file
     inquire(iolength=recl) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), &
          vkl(:,ik), xiou,xiuo
     inquire(iolength=reclfull) pwmat
     un=51
     un2=52
     open(unit=un,file='EMAT_Q00001.OUT',form='unformatted', &
          action='write',access='direct',recl=recl)
     write(un,rec=ik) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), vkl(:,ik), &
          xiou, xiuo
     close(un)
     ! write the matrix elements to direct-access file
!!$     write(50,rec=ik) pwmat
!!$     write(50,'(i8,3g18.10)') ik,vkl(:,ik)
!!$     do igq=1,ngq(iq)
!!$        do ist=1,nstsv
!!$           do jst=1,nstsv
!!$              write(50,'(4i8,3g18.10)') ik,igq,ist,jst,pwmat(igq,ist,jst), &
!!$                   abs(pwmat(igq,ist,jst))**2
!!$           end do
!!$        end do
!!$     end do
!!$     write(50,*)


     ! rotate matrix element if k-point set is reduced
     !
     ! M(G)(ka,q) = exp(-i(G+q)a^-1Ta) M(Ga^-1+G_a)(k,q) !

     if (nkpt.ne.nkptnr) then
        do j=1,nsymcrysstr(ik)
           iknr=ikstrmapiknr(j,ik)
           isym=scmapstr(j,ik)
           jsym=scimap(isym)
           lspl=lsplsymc(jsym)       
           ! rotation in Cartesian coordinates
           c(:,:)=symlatc(:,:,lspl)
           ! rotation in lattice coordinates
           s(:,:)=symlat(:,:,lspl)
           do igq=1,ngq(iq)
              t1=twopi*dot_product(vgql(:,igq,iq),matmul(s,vtlsymc(:,jsym)))
              t2=cos(t1); t3=-sin(t1)
              if (abs(t2).lt.epsrot) t2=0.d0
              if (abs(t3).lt.epsrot) t3=0.d0
              zt1=cmplx(t2,t3,8)
              ig=igqig(igq,iq)
              vg(:)=ivg(:,ig)
              vg=matmul(vg,s)
              ig=ivgigq(vg(1),vg(2),vg(3),iq)

              write(80,'(4i6,3x,3i5,3x,i6)') iknr,ik,isym,igq,vg,ig

              ! write to ASCII file
              do ist=1,nstsv
                 do jst=1,nstsv
                    zt2=zt1*pwmat(ig,ist,jst)
                    ! matrix elements for full k-point set
                    pwmatf(igq,ist,jst)=zt2
                    write(90,'(7i5,5g18.10)') iknr,ik,isym,igq,ist,jst,ig,zt2, &
                         abs(zt2)**2, zt1
                    
                    if ((ist.le.nstval).and.(jst.gt.nstval)) then
                       yiou(ist,jst-nstval,igq)=zt2
                       write(300 + iq,'(a,4i6,3g18.10)') 'ik,igq,i1,i2', &
                            iknr,igq,ist,jst-nstval,zt2,abs(zt2)**2
                    end if
                    if ((ist.gt.nstval).and.(jst.le.nstval)) then
                       yiuo(ist-nstval,jst,igq)=zt2
                       write(400 + iq,'(a,4i6,3g18.10)') 'ik,igq,i1,i2', &
                            iknr,igq,ist-nstval,jst,zt2,abs(zt2)**2
                    end if
                 end do
              end do
              ! end loop over G+q vectors
           end do

           ! write v-c and c-v matrix elements
           open(unit=un,file='EMAT_NR_Q00001.OUT',form='unformatted', &
                action='write',access='direct',recl=recl)
           write(un,rec=iknr) nstval, nstcon, nkptnr, ngq(iq), vql(:,iq), &
                vklnr(:,iknr), yiou, yiuo
           close(un)

           ! write full matrix elements
           open(unit=un,file='EMAT_FULL_NR_Q00001.OUT',form='unformatted', &
                action='write',access='direct',recl=reclfull)
           write(un,rec=iknr) pwmatf
           close(un)

           ! end loop over elements of star
        end do
     end if

     ! end loop over k
  end do
  close(50)

  ! read and write matrix elements of non-reduced k-points to ASCII file
  if (nkpt.ne.nkptnr) then
     open(unit=un,file='EMAT_FULL_NR_Q00001.OUT',form='unformatted', &
          action='read',access='direct',recl=reclfull)
     open(unit=un2,file='PWMAT_NR_ASC.OUT',form='formatted', &
          action='write',status='replace')
     do iknr=1,nkptnr
        read(un,rec=iknr) pwmat
        do igq=1,ngq(iq)
           do ist=1,nstsv
              do jst=1,nstsv
                 write(un2,'(4i5,3g18.10)') iknr,igq,ist,jst, &
                      pwmat(igq,ist,jst),abs(pwmat(igq,ist,jst))**2
              end do
           end do
        end do
     end do
     close(un2)
     close(un)
  end if

  write(*,*)
  write(*,'("Info(writepwmat):")')
  write(*,'(" matrix elements of the plane wave written to file PWMAT.OUT")')
  write(*,*)
  deallocate(apwalmk,evecfvk,evecsvk,apwalmkp,evecfvkp,evecsvkp,pwmat)
end subroutine writepwmat
!EOC
