
module m_dfqoscwg
  implicit none
contains

  subroutine dfqoscwg(sw,iq,ik,pou,puo,xou,xuo,you,yuo)
    use modmain
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: sw,iq,ik
    complex(8), intent(in) :: pou(3),puo(3),xou(:),xuo(:)
    complex(8), intent(out) :: you(:),yuo(:)
    ! local variables
    complex(8) :: pout,puot
    real(8) :: s1(3,3)
    integer :: oc,i

    integer :: j,n,isym,jsym,ig,igq,lspl,vg(3),s(3,3)
    real(8), parameter :: epsrot=1.d-12
    real(8) :: t1,t2,t3
    complex(8) :: zt1

    !//////// TO BE DONE /////////////////

    if (symwings) then
       s1(:,:)=0.5d0*(symdfq0(:,:)+transpose(symdfq0))
       pout=zzero
       puot=zzero
       do i=1,3
          pout=pout+s1(i,i)*pou(i)
          puot=pout+s1(i,i)*puo(i)
       end do
    else
       oc=optcomp(1,1)
       pout=pou(oc)
       puot=puo(oc)
    end if

!!$    if (sw.eq.1) then
!!$       ! first wing (G=0,G/=0)
!!$       you(:)=pout*conjg(xou(:))
!!$       yuo(:)=puot*conjg(xuo(:))
!!$    else
!!$       ! second wing (G/=0,G=0)
!!$       you(:)=xou(:)*conjg(pout)
!!$       yuo(:)=xuo(:)*conjg(puot)
!!$    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    if (nkpt.ne.nkptnr) then
if(.true.)then
       n=ngq(iq)
       you(:)=zzero
       yuo(:)=zzero
       do igq=2,ngq(iq)
          do j=1,nsymcrysstr(ik)
             isym=scmapstr(j,ik)
             jsym=scimap(isym)
             lspl=lsplsymc(jsym)       
             ! rotation in lattice coordinates
             s(:,:)=symlat(:,:,lspl)
             ! let c=a^(-1) and note that in exciting the symmetry operations
             ! are defined as y = as a (x + T_a)
             ! (G+q).T_c
             t1=twopi*dot_product(-vgql(:,igq,iq), &
                  matmul(s,vtlsymc(:,jsym)))
             ! exp(-i(G+q)T_c)
             t2=cos(t1); t3=-sin(t1)
             if (abs(t2).lt.epsrot) t2=0.d0
             if (abs(t3).lt.epsrot) t3=0.d0
             zt1=cmplx(t2,t3,8)
             ! G-vector of G+q
             ig=igqig(igq,iq)
             vg(:)=ivg(:,ig)
             ! G*c
             vg=matmul(vg,s)
             ! index of G+q-vector of new G-vector
             ig=ivgigq(vg(1),vg(2),vg(3),iq)
write(82,'(4i6,2g12.4)') ik,igq,j,ig,zt1
             ! summation over star of k-point
             you(igq-1)=you(igq-1)+zt1*pout*conjg(xou(ig-1))
             yuo(igq-1)=yuo(igq-1)+zt1*puot*conjg(xuo(ig-1))
          end do
       end do
    else
       ! first wing (G=0,G/=0)
       you(:)=pout*conjg(xou(:))
       yuo(:)=puot*conjg(xuo(:))
    end if
    if (sw.ne.1) then
       you=conjg(you)
       yuo=conjg(yuo)
    end if

  end subroutine dfqoscwg

end module m_dfqoscwg
