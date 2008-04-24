
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqkgmt(iq,ik,igq)
  use modmain
  use modxs
  use m_zaxpyc
  use m_tdzoutpr
  use m_tdzoutpr3
  implicit none
  ! arguments
  integer, intent(in) :: iq,ik,igq
  ! local variables
  character(*), parameter :: thisnam='ematqkgmt'
  integer :: is,ia,ias,l1,m1,lm1,l3,m3,lm3,io,io1,io2,ilo,ilo1,ilo2
  integer :: lmax1,lmax3, i,j,ikt
  complex(8), allocatable :: zv(:)
  allocate(zv(nstsv))
  ikt=ik
  lmax1=lmaxapwtd
  lmax3=lmax1
  xih(:,:) = zzero
  xiuhloa(:,:) = zzero
  xiohalo(:,:) = zzero
  xiou(:,:,igq)=zzero
  ! loop over species and atoms
  do is = 1, nspecies
     do ia = 1, natoms(is)
        ias = idxas(ia,is)
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
                               apwdlm(1,io2,lm3,ias),1,zv,1)
                       end do
                    end do ! m3
                 end do ! l3
                 call tdzoutpr(nst1,nst2, &
                      fourpi*conjg(sfacgq(igq,ias,iq)), &
                      apwdlm0(istlo1:isthi1,io1,lm1,ias),zv(istlo2:isthi2), &
                      xiou(:,:,igq))
                 ! end loop over (l',m',p')
              end do! io1
           end do ! m1
        end do ! l1
        call cpu_time(cmt1)


!!$        !--------------------------------------!
!!$        !     local-orbital-APW contribution   !
!!$        !--------------------------------------!
!!$        ! loop over local orbitals
!!$        do ilo=1,nlorb(is)
!!$           l1=lorbl(ilo,is)
!!$           do m1=-l1,l1
!!$              lm1=idxlm(l1,m1)
!!$              i=idxlo(lm1,ilo,ias)
!!$              zv(:)=zzero
!!$              ! loop over (l'',m'',p'')
!!$              do l3=0,lmax3
!!$                 do m3=-l3,l3
!!$                    lm3=idxlm(l3,m3)
!!$                    do io=1,apword(l3,is)
!!$                       call zaxpy(nstsv, &
!!$                            intrgloa(lm1,ilo,lm3,io,ias), &
!!$                            apwdlm(1,io,lm3,ias),1,zv,1)
!!$                    end do ! io
!!$                 end do ! m3
!!$              end do ! l3
!!$              xiuhloa(i,:)=xiuhloa(i,:)+fourpi*conjg(sfacgq(igq,ias,iq))* &
!!$                   zv(istlo2:isthi2)
!!$           end do ! m1
!!$        end do ! ilo

        !--------------------------------------!
        !     local-orbital-APW contribution   ! new
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
                            intrgloa(lm1,ilo,lm3,io,ias), &
                            apwdlm(1,io,lm3,ias),1,zv,1)
                    end do ! io
                 end do ! m3
              end do ! l3
              call tdzoutpr(nst1,nst2, &
                   fourpi*conjg(sfacgq(igq,ias,iq)), &
                   lodlm0(istlo1:isthi1,ilo,lm1,ias),zv(istlo2:isthi2), &
                   xiou(:,:,igq))
           end do ! m1
        end do ! ilo

        call cpu_time(cmt2)


!!$        !--------------------------------------!
!!$        !     APW-local-orbital contribution   !
!!$        !--------------------------------------!
!!$        ! loop over local orbitals
!!$        do ilo=1,nlorb(is)
!!$           l1=lorbl(ilo,is)
!!$           do m1=-l1,l1
!!$              lm1=idxlm(l1,m1)
!!$              j=idxlo(lm1,ilo,ias)
!!$              zv(:)=zzero
!!$              ! loop over (l'',m'',p'')
!!$              do l3=0,lmax3
!!$                 do m3=-l3,l3
!!$                    lm3=idxlm(l3,m3)
!!$                    do io=1,apword(l3,is)
!!$                       call zaxpyc(nstsv, &
!!$                            intrgalo(lm1,ilo,lm3,io,ias), &
!!$                            apwdlm0(:,io,lm3,ias),1,zv,1)
!!$                    end do ! io
!!$                 end do ! m3
!!$              end do ! l3
!!$              xiohalo(:,j)=xiohalo(:,j)+fourpi*conjg(sfacgq(igq,ias,iq))* &
!!$                   zv(istlo1:isthi1)
!!$           end do ! m1
!!$        end do ! ilo

!!$        !--------------------------------------!
!!$        !     APW-local-orbital contribution   ! new fail
!!$        !--------------------------------------!
!!$        ! loop over local orbitals
!!$        do ilo=1,nlorb(is)
!!$           l1=lorbl(ilo,is)
!!$           do m1=-l1,l1
!!$              lm1=idxlm(l1,m1)
!!$              zv(:)=zzero
!!$              ! loop over (l'',m'',p'')
!!$              do l3=0,lmax3
!!$                 do m3=-l3,l3
!!$                    lm3=idxlm(l3,m3)
!!$                    do io=1,apword(l3,is)
!!$                       call zaxpyc(nstsv, &
!!$                            intrgalo(lm1,ilo,lm3,io,ias), &
!!$                            apwdlm0(:,io,lm3,ias),1,zv,1)
!!$                    end do ! io
!!$                 end do ! m3
!!$              end do ! l3
!!$              call tdzoutpr3(nst1,nst2, &
!!$                   fourpi*conjg(sfacgq(igq,ias,iq)), &
!!$                   lodlm(istlo1:isthi1,ilo,lm1,ias),zv(istlo2:isthi2), &
!!$                   xiou(:,:,igq))
!!$           end do ! m1
!!$        end do ! ilo

        !--------------------------------------!
        !     APW-local-orbital contribution   ! new
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
                            intrgalo(lm1,ilo,lm3,io,ias), &
                            lodlm(1,ilo,lm1,ias),1,zv,1)
                    end do ! m1
                 end do ! ilo
                 call tdzoutpr(nst1,nst2, &
                      fourpi*conjg(sfacgq(igq,ias,iq)), &
                      apwdlm0(istlo1:isthi1,io,lm3,ias),zv(istlo2:isthi2), &
                      xiou(:,:,igq))
              end do ! io
           end do ! m3
        end do ! l3
        call cpu_time(cmt3)

!!$        !------------------------------------------------!
!!$        !     local-orbital-local-orbital contribution   !
!!$        !------------------------------------------------!
!!$        do ilo1=1,nlorb(is)
!!$           l1=lorbl(ilo1,is)
!!$           do m1=-l1,l1
!!$              lm1=idxlm(l1,m1)
!!$              i=idxlo(lm1,ilo1,ias)
!!$              do ilo2=1,nlorb(is)
!!$                 l3=lorbl(ilo2,is)
!!$                 do m3=-l3,l3
!!$                    lm3=idxlm(l3,m3)
!!$                    j=idxlo(lm3,ilo2,ias)
!!$                    xih(i,j) = xih(i,j) + &
!!$                         fourpi*conjg(sfacgq(igq,ias,iq))* &
!!$                         intrglolo(lm1,ilo1,lm3,ilo2,ias)
!!$                 end do ! m3
!!$              end do ! ilo2
!!$           end do ! m1
!!$        end do ! ilo1


        !------------------------------------------------!
        !     local-orbital-local-orbital contribution   ! new
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
                         intrglolo(lm1,ilo1,lm3,ilo2,ias), &
                         lodlm(1,ilo2,lm3,ias),1,zv,1)
                 end do ! m3
              end do ! ilo2
              call tdzoutpr(nst1,nst2, &
                   fourpi*conjg(sfacgq(igq,ias,iq)), &
                   lodlm0(istlo1:isthi1,ilo1,lm1,ias),zv(istlo2:isthi2), &
                   xiou(:,:,igq))
           end do ! m1
        end do ! ilo1


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
