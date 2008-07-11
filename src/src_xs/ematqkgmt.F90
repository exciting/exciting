
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqkgmt(iq,ik,igq)
  use modmain
  use modxs
  use m_zaxpyc
  use m_xszoutpr
  use m_xszoutpr3
  implicit none
  ! arguments
  integer, intent(in) :: iq,ik,igq
  ! local variables
  character(*), parameter :: thisnam='ematqkgmt'
  integer :: is,ia,ias,l1,m1,lm1,l3,m3,lm3,io,io1,io2,ilo,ilo1,ilo2
  integer :: lmax1,lmax3,ikt,i,j
  complex(8), allocatable :: zv(:)
  allocate(zv(nstsv))
  ikt=ik
  lmax1=lmaxapwwf
  lmax3=lmax1
  xih(:,:) = zzero
  xiuhloa(:,:) = zzero
  xiohalo(:,:) = zzero
  xiou(:,:,igq)=zzero
  ! loop over species and atoms
  do is=1,nspecies
     do ia=1,natoms(is)
        ias=idxas(ia,is)
        call cpu_time(cmt0)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
        ! loop over (l',m',p')
        do l1=0,lmax1
           do m1=-l1,l1
              lm1=idxlm(l1,m1)
              do io1=1,apword(l1,is)
                 zv(:)=zzero
                 ! loop over (l'',m'',p'')
                 do l3=0,lmax3
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       do io2=1,apword(l3,is)
                          call zaxpy(nstsv, &
                               intrgaa(lm1,io1,lm3,io2,ias), &
                               apwcmt(1,io2,lm3,ias),1,zv,1)
                       end do
                    end do ! m3
                 end do ! l3
                 call xszoutpr(nst1,nst2, &
                      fourpi*conjg(sfacgq(igq,ias,iq)), &
                      apwcmt0(istlo1:isthi1,io1,lm1,ias),zv(istlo2:isthi2), &
                      xiou(:,:,igq))
                 ! end loop over (l',m',p')
              end do! io1
           end do ! m1
        end do ! l1
        call cpu_time(cmt1)
        if (fastemat) then
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
           ! loop over local orbitals
           do ilo=1,nlorb(is)
              l1=lorbl(ilo,is)
              do m1=-l1,l1
                 lm1=idxlm(l1,m1)
                 zv(:)=zzero
                 ! loop over (l'',m'',p'')
                 do l3=0,lmax3
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       do io=1,apword(l3,is)
                          call zaxpy(nstsv, &
                               intrgloa(m1,ilo,lm3,io,ias), &
                               apwcmt(1,io,lm3,ias),1,zv,1)
                       end do ! io
                    end do ! m3
                 end do ! l3
                 call xszoutpr(nst1,nst2, &
                      fourpi*conjg(sfacgq(igq,ias,iq)), &
                      locmt0(istlo1:isthi1,ilo,m1,ias),zv(istlo2:isthi2), &
                      xiou(:,:,igq))
              end do ! m1
           end do ! ilo
           call cpu_time(cmt2)
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
           ! loop over (l'',m'',p'')
           do l3=0,lmax3
              do m3=-l3,l3
                 lm3=idxlm(l3,m3)
                 do io=1,apword(l3,is)
                    zv(:)=zzero
                    ! loop over local orbitals
                    do ilo=1,nlorb(is)
                       l1=lorbl(ilo,is)
                       do m1=-l1,l1
                          lm1=idxlm(l1,m1)
                          call zaxpy(nstsv, &
                               intrgalo(m1,ilo,lm3,io,ias), &
                               locmt(1,ilo,m1,ias),1,zv,1)
                       end do ! m1
                    end do ! ilo
                    call xszoutpr(nst1,nst2, &
                         fourpi*conjg(sfacgq(igq,ias,iq)), &
                         apwcmt0(istlo1:isthi1,io,lm3,ias),zv(istlo2:isthi2), &
                         xiou(:,:,igq))
                 end do ! io
              end do ! m3
           end do ! l3
           call cpu_time(cmt3)
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
           do ilo1=1,nlorb(is)
              l1=lorbl(ilo1,is)
              do m1=-l1,l1
                 lm1=idxlm(l1,m1)
                 zv(:)=zzero
                 do ilo2=1,nlorb(is)
                    l3=lorbl(ilo2,is)
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       call zaxpy(nstsv, &
                            intrglolo(m1,ilo1,m3,ilo2,ias), &
                            locmt(1,ilo2,m3,ias),1,zv,1)
                    end do ! m3
                 end do ! ilo2
                 call xszoutpr(nst1,nst2, &
                      fourpi*conjg(sfacgq(igq,ias,iq)), &
                      locmt0(istlo1:isthi1,ilo1,m1,ias),zv(istlo2:isthi2), &
                      xiou(:,:,igq))
              end do ! m1
           end do ! ilo1
        else
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
           ! loop over local orbitals
           do ilo=1,nlorb(is)
              l1=lorbl(ilo,is)
              do m1=-l1,l1
                 lm1=idxlm(l1,m1)
                 i=idxlo(lm1,ilo,ias)
                 zv(:)=zzero
                 ! loop over (l'',m'',p'')
                 do l3=0,lmax3
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       do io=1,apword(l3,is)
                          call zaxpy(nstsv, &
                               intrgloa(m1,ilo,lm3,io,ias), &
                               apwcmt(1,io,lm3,ias),1,zv,1)
                       end do ! io
                    end do ! m3
                 end do ! l3
                 xiuhloa(i,:)=xiuhloa(i,:)+fourpi*conjg(sfacgq(igq,ias,iq))* &
                      zv(istlo2:isthi2)
              end do ! m1
           end do ! ilo
           call cpu_time(cmt2)
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
           ! loop over local orbitals
           do ilo=1,nlorb(is)
              l1=lorbl(ilo,is)
              do m1=-l1,l1
                 lm1=idxlm(l1,m1)
                 j=idxlo(lm1,ilo,ias)
                 zv(:)=zzero
                 ! loop over (l'',m'',p'')
                 do l3=0,lmax3
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       do io=1,apword(l3,is)
                          call zaxpyc(nstsv, &
                               intrgalo(m1,ilo,lm3,io,ias), &
                               apwcmt0(:,io,lm3,ias),1,zv,1)
                       end do ! io
                    end do ! m3
                 end do ! l3
                 xiohalo(:,j)=xiohalo(:,j)+fourpi*conjg(sfacgq(igq,ias,iq))* &
                      zv(istlo1:isthi1)
              end do ! m1
           end do ! ilo
           call cpu_time(cmt3)
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
           do ilo1=1,nlorb(is)
              l1=lorbl(ilo1,is)
              do m1=-l1,l1
                 lm1=idxlm(l1,m1)
                 i=idxlo(lm1,ilo1,ias)
                 do ilo2=1,nlorb(is)
                    l3=lorbl(ilo2,is)
                    do m3=-l3,l3
                       lm3=idxlm(l3,m3)
                       j=idxlo(lm3,ilo2,ias)
                       xih(i,j) = xih(i,j) + &
                            fourpi*conjg(sfacgq(igq,ias,iq))* &
                            intrglolo(m1,ilo1,m3,ilo2,ias)
                    end do ! m3
                 end do ! ilo2
              end do ! m1
           end do ! ilo1
           ! strategy of emat-calculation
        end if
        call cpu_time(cmt4)
        cpumtaa=cpumtaa+cmt1-cmt0
        cpumtloa=cpumtloa+cmt2-cmt1
        cpumtalo=cpumtalo+cmt3-cmt2
        cpumtlolo=cpumtlolo+cmt4-cmt3
        ! end loop over species and atoms
     end do ! ia
  end do ! is
  ! deallocate
  deallocate(zv)
end subroutine ematqkgmt
