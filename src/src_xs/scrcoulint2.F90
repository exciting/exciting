
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

integer function iplocnr(ivp,ngridp)
  ! locate p-point with index on grid following the convention that the
  ! first coordinate runs fastest.
  implicit none
  ! arguments
  integer, intent(in) :: ivp(3),ngridp(3)
  ! this should be consistent with ipmap of "genppts.f90"
  iplocnr = 1 + ivp(1) + ngridp(1)*ivp(2) + ngridp(1)*ngridp(2)*ivp(3)
end function iplocnr


!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////


subroutine getscreen(iqr,ngq,scrh,scrw,scrb)
  use m_genfilname
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iqr,ngq
  complex(8), intent(out) :: scrb(ngq,ngq),scrw(ngq,2,3),scrh(9)
  ! local variables
  character(256) :: fname
  real(8) :: rm(2,9)
  integer :: igq1,igq2,j,it1,it2,it3,un
  ! read in screening
  call genfilname(basename='SCREEN',iq=iqr,filnam=fname)
  call getunit(un)
  open(un,file=trim(fname),form='formatted',action='read',status='old')
  do igq1=1,ngq
     do igq2=1,ngq
        if (iqr.eq.1) then
           if ((igq1.eq.1).and.(igq2.eq.1)) then
              read(un,*) (it1,it2,it3,rm(1,j),rm(2,j),j=1,9)
              scrh(:)=cmplx(rm(1,:),rm(2,:),8)
           end if
           if ((igq1.eq.1).and.(igq2.ne.1)) then
              read(un,*) (it1,it2,it3,rm(1,j),rm(2,j),j=1,3)
              scrw(igq2,1,:)=cmplx(rm(1,:3),rm(2,:3),8)
           end if
           if ((igq1.ne.1).and.(igq2.eq.1)) then
              read(un,*) (it1,it2,it3,rm(1,j),rm(2,j),j=1,3)
              scrw(igq1,2,:)=cmplx(rm(1,:3),rm(2,:3),8)
           end if
           if ((igq1.ne.1).and.(igq2.ne.1)) read(un,*) it1,it2,it3,&
                rm(1,1),rm(2,1)
           scrb(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
        else
           read(un,*) it1,it2,it3,rm(1,1),rm(2,1)
           scrb(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
        end if
     end do
  end do
  close(un)
end subroutine getscreen


!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////


subroutine geniscreen(iqr,nmax,n,scri)
  use modxs
  use invert
  implicit none
  ! arguments
  integer, intent(in) :: iqr,n,nmax
  complex(8), intent(out) :: scri(nmax,nmax)
  ! local variables
  complex(8), allocatable :: tm(:,:),tmi(:,:)
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:)
  complex(8) :: scrnh0(3),scrnih0(3)
  integer :: oct,j
  integer, external :: octmap

  allocate(scrn(n,n),scrnw(n,2,3),scrnh(9),tm(n,n),tmi(n,n))
  ! read screening from file
  call getscreen(iqr,n,scrnh,scrnw,scrn)

  ! invert dielectric matrix for q-point going to zero in x-, y-, and
  ! z-direction for q=0
  if (iqr.eq.1) then
     do oct=1,3
        ! index for diagonal tensor component
        j=octmap(oct,oct)
        scrn(1,1)=scrnh(j)
        ! keep head of dielectric matrix for q=0
        scrnh0(oct)=scrn(1,1)
        if (n.gt.1) then
           scrn(1,2:n)=scrnw(2:n,1,oct)
           scrn(2:n,1)=scrnw(2:n,2,oct)
        end if
        ! store screening in temporary array
        tm=scrn

write(*,'(a,i5,2g18.10)') 'optcomp,    1/eps_00(q=0)      :', &
             oct,1.d0/scrnh0(oct)
write(*,'(a,i5,2g18.10)') 'optcomp,      eps_00(q=0)      :', &
             oct,scrnh0(oct)

        call zinvert_hermitian(scrherm,tm,tmi)
        scri(1:n,1:n)=tmi(:,:)
        ! keep head of inverse dielectric matrix for q=0
        scrnih0(oct)=scri(1,1)

write(*,'(a,i5,2g18.10)') 'optcomp,   eps^-1_00(q=0)      :', &
     oct,scrnih0(oct)
write(*,'(a,i5,2g18.10)') 'optcomp, 1/eps^-1_00(q=0)      :', &
     oct,1.d0/scrnih0(oct)

     end do
     ! symmetrize head of inverse dielectric matrix wrt. the directions
     ! in which q goes to zero
     call symsci0(bsediagsym,scrnh0,scrnih0,scri(1,1))

write(*,'(a,i5,2g18.10)') 'optcomp, symm.   eps^-1_00(q=0):', &
     iqr,scri(1,1)
write(*,'(a,i5,2g18.10)') 'optcomp, symm. 1/eps^-1_00(q=0):', &
          iqr,1.d0/scri(1,1)

  else
     tm(:,:)=scrn(:,:)

write(*,'(a,i5,2g18.10)') 'iq,    1/eps_00(q)             :', &
     iqr,1.d0/scrn(1,1)
write(*,'(a,i5,2g18.10)') 'iq,      eps_00(q)             :', &

     iqr,scrn(1,1)
     call zinvert_hermitian(scrherm,tm,tmi)
     scri(1:n,1:n)=tmi(:,:)
  end if

write(*,'(a,i5,2g18.10)') 'diel. matr.: iq,   eps^-1_00(q):', &
     iqr,scri(1,1)
write(*,'(a,i5,2g18.10)') 'diel. matr.: iq, 1/eps^-1_00(q):', &
     iqr,1.d0/scri(1,1)
write(*,*)

  deallocate(scrn,scrnw,scrnh,tm,tmi)
end subroutine geniscreen


!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////


subroutine findsymeqiv(vpl,vplr,nsc,sc,ivgsc)
  use modmain
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3),vplr(3)
  integer, intent(out) :: nsc,sc(maxsymcrys),ivgsc(3,maxsymcrys)
  ! local variables
  integer :: isym,lspl,iv(3)
  real(8) :: s(3,3),v1(3),t1
  real(8), external :: r3taxi
  ! symmetries that transform non-reduced q-point to reduced one, namely
  ! q1 = s^-1 * q + G_s. Here, q1 is vpl, q is vplr.
  nsc=0
  do isym=1,nsymcrys
     lspl=lsplsymc(isym)
     s(:,:)=dble(symlat(:,:,lspl))
     call r3mtv(s,vplr,v1)
     call r3frac(epslat,v1,iv)
     t1=r3taxi(vpl,v1)
     if (t1.lt.epslat) then
        nsc=nsc+1
        sc(nsc)=isym
        ivgsc(:,nsc)=-iv(:)
     end if
  end do
  if (nsc.eq.0) then
     write(*,*)
     write(*,'("Error(findsymeqiv): p-points are not equivalent by symmetry")')
     write(*,'(" vpl  :",3g18.10)') vpl
     write(*,'(" vplr :",3g18.10)') vplr
     write(*,*)
     call terminate
  end if
end subroutine findsymeqiv


!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////


subroutine findgqmap(iq,iqr,nsc,sc,ivgsc,nmax,n,isc,isci,ivgu,igqmap)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iq,iqr,nsc,sc(maxsymcrys),ivgsc(3,maxsymcrys),nmax,n
  integer, intent(out) :: isc,isci,ivgu(3),igqmap(nmax)
  ! local variables
  real(8) :: vqr(3),v2(3),t1
  integer :: iqrnr,j,isym,isymi,lspl,lspli,iv(3),ivg1(3),igq1
  integer, external :: iplocnr
  ! find map from G-vectors to rotated G-vectors
  iqrnr=iplocnr(ivqr(1,iqr),ngridq)
  vqr(:)=vqlr(:,iqr)
  do j=1,nsc
     isym=sc(j)
     lspl=lsplsymc(isym)
     isymi=scimap(isym)
     lspli=lsplsymc(isymi)
     do igq1=1,n
        ivg1(:)=ivg(:,igqig(igq1,iq))
        ! G1 = si^-1 * ( G + G_s ) , where si is the inverse of s
        iv=matmul(transpose(symlat(:,:,lspli)),ivg1+ivgsc(:,j))
        ! |G1 + q|
        v2=matmul(bvec,iv+vqr)
        t1=sqrt(sum(v2**2))
!!$write(*,'(a,5i6,3g18.10,3x,3g18.10)') 'reduce:',iq,j,isym,&
!!$     igq1,ivgigq(iv(1),iv(2),iv(3),iqrnr),vgql(:,igq1,iq)- &
!!$     matmul(transpose(symlat(:,:,lspl)),vqr+iv)
        if ((n.gt.1).and.(t1.gt.gqmax)) then
           write(*,*) '*** need one more symmetry operation'
           goto 10
        end if
        ! locate G1 + q in G+q-vector set
        igqmap(igq1)=ivgigq(iv(1),iv(2),iv(3),iqrnr)
        if (igqmap(igq1).le.0) then
           write(*,*)
           write(*,'("Error(findgqmap): failed to map rotated G-vector")')
           write(*,'(" non-reduced q-point                    :",i8)') iq
           write(*,'(" reduced q-point                        :",i8)') iqr
           write(*,'(" reduced q-point in non-reduced set     :",i8)') iqrnr
           write(*,'(" G+q-vector index (non-reduced q-point) :",i8)') igq1
           write(*,'(" rotated G-vector                       :",3i8)') iv
           write(*,*)
           call terminate
        end if
        ! end loop over G
     end do
     ! store G1 vector
     ivgu(:)=ivgsc(:,j)
     isc=isym
     isci=isymi
     goto 20
10   continue
     ! end loop over symmetry operations
  end do
  write(*,*)
  write(*,'("Error(findgqmap): failed to reduce q-point: ",i8)') iq
  write(*,*)
  call terminate
20 continue
end subroutine findgqmap

!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////

subroutine genphasedm(iq,jsym,nmax,n,phfdm)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iq,jsym,nmax,n
  complex(8), intent(out) :: phfdm(nmax,nmax)
  ! local variables
  real(8), parameter :: epsortho=1.d-12
  real(8) :: vtl(3),t1,t2,t3,s(3,3),si(3,3)
  integer :: igq1,igq2,ivg1(3),ivg2(3),iv(3),jsymi
  jsymi=scimap(jsym)
  s(:,:)=dble(symlat(:,:,lsplsymc(jsym)))
  si(:,:)=dble(symlat(:,:,lsplsymc(jsymi)))

  do igq1=1,n
     ivg1(:)=ivg(:,igqig(igq1,iq))
     do igq2=igq1,n
        ! G-vector difference
        ivg2(:)=ivg(:,igqig(igq2,iq))-ivg1(:)

        !*** try change in sign of G-vector difference **********************
        ivg2=-ivg2

        ! translation vector s^-1*vtl(s^-1)
        !                 vtl=matmul(transpose(s),vtlsymc(:,jsymi))
        ! we multiply in real space
        vtl=matmul(si,vtlsymc(:,jsymi))
        call r3frac(epslat,vtl,iv)
        t1=twopi*dot_product(dble(ivg2),vtl)
        t2=cos(t1)
        t3=sin(t1)
        if (abs(t2).lt.epsortho) t2=0.d0
        if (abs(t3).lt.epsortho) t3=0.d0
        ! phase factor for dielectric matrix (due to translations)
        phfdm(igq1,igq2)=cmplx(t2,t3,8)
        phfdm(igq2,igq1)=conjg(phfdm(igq1,igq2))

!write(40,'(a,i5,2x,2i5,2x,2i5,2g18.10)') 'q,g,gp,isym,isymi,phf',iq,igq1,igq2,jsym,jsymi,phfdm(igq1,igq2)

        ! end loop over (G,Gp)-vectors
     end do
  end do
end subroutine genphasedm




