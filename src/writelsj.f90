
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writelsj
use modmain
implicit none
! local variables
integer ik,ispn,jspn,i,j,kst
integer is,ia,ias,ist,irc,n
real(8) xl(3),xs(3),xj(3),wo
complex(8) zt1
! allocatable arrays
logical, allocatable :: done(:,:)
real(8), allocatable :: xlt(:,:),xst(:,:),xjt(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wfmt3(:,:,:,:)
complex(8), allocatable :: wfmt4(:,:,:,:)
complex(8), allocatable :: zlflm(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
! initialise universal variables
call init0
call init1
allocate(done(nstfv,nspnfv))
allocate(xlt(3,natmtot),xst(3,natmtot),xjt(3,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt1(lmmaxapw,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxapw,nrcmtmax,2))
allocate(wfmt3(lmmaxapw,nrcmtmax,2,3))
allocate(wfmt4(lmmaxapw,nrcmtmax,2,3))
allocate(zlflm(lmmaxapw,3))
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
if (task.eq.15) then
! get the occupancies from file
  do ik=1,nkpt
    call getoccsv(vkl(1,ik),occsv(1,ik))
  end do
else
! flag the occupancies for the k-points and states in kstlist
  occsv(:,:)=0.d0
  do kst=1,nkstlist
    ik=kstlist(1,kst)
    j=kstlist(2,kst)
    if ((ik.le.0).or.(ik.gt.nkpt)) then
      write(*,*)
      write(*,'("Error(writelsj): k-point out of range : ",I8)') ik
      write(*,*)
      stop
    end if
    if ((j.le.0).or.(j.gt.nstsv)) then
      write(*,*)
      write(*,'("Error(writelsj): state out of range : ",I8)') j
      write(*,*)
      stop
    end if
    occsv(j,ik)=1.d0
  end do
end if
if (task.eq.15) then
  open(50,file='LSJ.OUT',action='WRITE',form='FORMATTED')
else
  open(50,file='LSJ_KST.OUT',action='WRITE',form='FORMATTED')
end if
write(50,*)
write(50,'("Expectation values are computed only over the muffin-tin")')
! zero the total expection values
xlt(:,:)=0.d0
xst(:,:)=0.d0
xjt(:,:)=0.d0
! begin loop over k-point and state list
do ik=1,nkpt
! get the eigenvectors from file
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
     sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
  end do
! begin loop over atoms and species
  do is=1,nspecies
    n=lmmaxapw*nrcmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      done(:,:)=.false.
! begin loop over second-variational states
      do j=1,nstsv
        wo=wkpt(ik)*occsv(j,ik)
        if (abs(wo).gt.epsocc) then
          if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
            wfmt2(:,:,:)=0.d0
            i=0
            do ispn=1,nspinor
              if (spinsprl) then
                jspn=ispn
              else
                jspn=1
              end if
              do ist=1,nstfv
                i=i+1
                zt1=evecsv(i,j)
                if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                  if (.not.done(ist,jspn)) then
                    call wavefmt(lradstp,lmaxapw,is,ia,ngk(ik,jspn), &
                     apwalm(1,1,1,1,jspn),evecfv(1,ist,jspn),lmmaxapw, &
                     wfmt1(1,1,ist,jspn))
                    done(ist,jspn)=.true.
                  end if
! add to spinor wavefunction
                  call zaxpy(n,zt1,wfmt1(1,1,ist,jspn),1,wfmt2(1,1,ispn),1)
                end if
              end do
            end do
          else
! spin-unpolarised wavefunction
            call wavefmt(lradstp,lmaxapw,is,ia,ngk(ik,1),apwalm,evecfv(1,j,1), &
             lmmaxapw,wfmt2)
          end if
! apply the L operator
          do ispn=1,nspinor
            do irc=1,nrcmt(is)
              call lopzflm(lmaxapw,wfmt2(1,irc,ispn),lmmaxapw,zlflm)
              do i=1,3
                wfmt3(:,irc,ispn,i)=zlflm(:,i)
              end do
            end do
          end do
! compute the expectation value of L
          xl(:)=0.d0
          do i=1,3
            do ispn=1,nspinor
              zt1=zfmtinp(.true.,lmaxapw,nrcmt(is),rcmt(1,is),lmmaxapw, &
               wfmt2(1,1,ispn),wfmt3(1,1,ispn,i))
              xl(i)=xl(i)+dble(zt1)
            end do
          end do
! compute the expecations value of S
          if (spinpol) then
! apply the S operator
            do irc=1,nrcmt(is)
! S_x
              wfmt4(:,irc,1,1)=0.5d0*wfmt2(:,irc,2)
              wfmt4(:,irc,2,1)=0.5d0*wfmt2(:,irc,1)
! S_y
              wfmt4(:,irc,1,2)=-0.5d0*zi*wfmt2(:,irc,2)
              wfmt4(:,irc,2,2)=0.5d0*zi*wfmt2(:,irc,1)
! S_z
              wfmt4(:,irc,1,3)=0.5d0*wfmt2(:,irc,1)
              wfmt4(:,irc,2,3)=-0.5d0*wfmt2(:,irc,2)
            end do
            xs(:)=0.d0
            do i=1,3
              do ispn=1,nspinor
                zt1=zfmtinp(.true.,lmaxapw,nrcmt(is),rcmt(1,is),lmmaxapw, &
                 wfmt2(1,1,ispn),wfmt4(1,1,ispn,i))
                xs(i)=xs(i)+dble(zt1)
              end do
            end do
          else
! spin-unpolarised case
            xs(:)=0.d0
          end if
! expectation value of J
          xj(:)=xl(:)+xs(:)
! add to the total
          xlt(:,ias)=xlt(:,ias)+wo*xl(:)
          xst(:,ias)=xst(:,ias)+wo*xs(:)
          xjt(:,ias)=xjt(:,ias)+wo*xj(:)
! write the particular L, S and J to file
          if (task.eq.16) then
            write(50,*)
            write(50,'("k-point : ",I6,3G18.10)') ik,vkl(:,ik)
            write(50,'("state : ",I6)') j
            write(50,'("species : ",I4," (",A,"), atom : ",I4)') is, &
             trim(spsymb(is)),ia
            write(50,'(" L : ",3G18.10)') xl
            write(50,'(" S : ",3G18.10)') xs
            write(50,'(" J : ",3G18.10)') xj
          end if
        end if
! end loop over states
      end do
! end loops over atoms and species
    end do
  end do
! end loop over k-points
end do
! write the total L, S and J to file if required
if (task.eq.15) then
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
       trim(spsymb(is)),ia
      write(50,'(" L : ",3G18.10)') xlt(:,ias)
      write(50,'(" S : ",3G18.10)') xst(:,ias)
      write(50,'(" J : ",3G18.10)') xjt(:,ias)
    end do
  end do
end if
close(50)
write(*,*)
write(*,'("Info(writelsj):")')
if (task.eq.15) then
  write(*,'(" total L, S and J expectation values written to LSJ.OUT")')
else
  write(*,'(" L, S and J expectation values for each k-point and state")')
  write(*,'("  in kstlist written to LSJ_KST.OUT")')
end if
write(*,*)
deallocate(done,xlt,xst,xjt)
deallocate(apwalm,evecfv,evecsv)
deallocate(wfmt1,wfmt2,wfmt3,wfmt4,zlflm)
return
end subroutine

