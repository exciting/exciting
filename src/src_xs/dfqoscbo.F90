
module m_dfqoscbo
  implicit none
contains

  subroutine dfqoscbo(iq,ik,n,xou,xuo,you,yuo)
    use modmain
    use modxs
    use m_tdzoutpr2
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,n
    complex(8), intent(in) :: xou(:),xuo(:)
    complex(8), intent(out) :: you(:,:),yuo(:,:)
    ! local variables
    complex(8) :: zt

    integer :: a,b,j,isym,jsym,ig,igp,igq,igqp,lspl,vg(3),vgp(3),s(3,3)
    real(8), parameter :: epsrot=1.d-12
    real(8) :: t1,t2,t3
    complex(8) :: zt1

!!$    if (nkpt.ne.nkptnr) then
if(.true.) then
       a=1
       b=0
       if (n.lt.ngq(iq)) then
          a=2
          b=1
       end if
       you(:,:)=zzero
       yuo(:,:)=zzero
       do igq=a,ngq(iq)
          do igqp=a,ngq(iq)
             do j=1,nsymcrysstr(ik)
                isym=scmapstr(j,ik)
                jsym=scimap(isym)
                lspl=lsplsymc(jsym)       
                ! rotation in lattice coordinates
                s(:,:)=symlat(:,:,lspl)
                ! let c=a^(-1) and note that in exciting the symmetry operations
                ! are defined as y = as a (x + T_a)
                ! (G+q).T_c
                t1=twopi*dot_product(vgql(:,igq,iq)-vgql(:,igqp,iq), &
                     matmul(s,vtlsymc(:,jsym)))
                ! exp(-i(G+q)T_c)
                t2=cos(t1); t3=-sin(t1)
                if (abs(t2).lt.epsrot) t2=0.d0
                if (abs(t3).lt.epsrot) t3=0.d0
                zt1=cmplx(t2,t3,8)
                ! *** G-vector of G+q
                ig=igqig(igq,iq)
                vg(:)=ivg(:,ig)
                ! G*c
                vg=matmul(vg,s)
                ! index of G+q-vector of new G-vector
                ig=ivgigq(vg(1),vg(2),vg(3),iq)
                ! *** Gp-vector of Gp+q
                igp=igqig(igqp,iq)
                vgp(:)=ivg(:,igp)
                ! Gp*c
                vgp=matmul(vgp,s)
                ! index of Gp+q-vector of new Gp-vector
                igp=ivgigq(vgp(1),vgp(2),vgp(3),iq)
write(81,'(6i6,2g12.4)') ik,igq,igqp,j,ig,igp,zt1
                ! summation over star of k-point
                you(igq-b,igqp-b)=you(igq-b,igqp-b)+zt1*xou(ig-b)* &
                     conjg(xou(igp-b))
                yuo(igq-b,igqp-b)=yuo(igq-b,igqp-b)+zt1*xuo(ig-b)* &
                     conjg(xuo(igp-b))
             end do
          end do
       end do
    else
       zt=(1.d0,0.d0)
       call tdzoutpr2(n,n,zt,xou,xou,you)
       call tdzoutpr2(n,n,zt,xuo,xuo,yuo)
    end if
    
  end subroutine dfqoscbo

end module m_dfqoscbo
