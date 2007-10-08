
module m_tdlinopt
  implicit none
contains

  subroutine tdlinopt(iq)
    use modmain
    use modtddft
    use modtetra
    use modpar
    use m_findexciton
    use m_genwgrid
    use m_pade
    use m_genloss
    use m_gensigma
    use m_gensumrls
    use m_gensymdf
    use m_writeeps
    use m_writeloss
    use m_writesigma
    use m_writesumrls
    use m_writeexciton
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'tdlinopt'
    character(256) :: filnam,filnam2,str,fnexciton2
    complex(8),allocatable :: mdf(:), mdf1(:),w(:),wr(:),sigma(:)
    real(8),allocatable :: wplot(:),loss(:)
    real(8),allocatable :: eps1(:),eps2(:),cf(:,:)
    real(8) :: t1,wint(2),sumrls(3),brd
    complex(8) :: zt1
    integer :: n,m,recl,i,j,iw,wi,wf,nwdfp,iwt,nc,oct,optcompt(3)
    logical :: tq0

    tq0 = tq1gamma.and.(iq.eq.1)
    ! number of components (3 for q=0)
    nc=1
    if (tq0) nc=3

    ! file extension for q-point
    write(filext,'("_Q",i5.5,".OUT")') iq

    ! limits for w-points
    wi=wpari
    wf=wparf
    nwdfp=wparf-wpari+1

    ! matrix size for local field effects
    n=ngq(iq)
    allocate(mdf1(nwdf),w(nwdf),wr(nwdos),wplot(nwdos),mdf(nwdos), &
         loss(nwdos),sigma(nwdos),cf(3,nwdos))
    allocate(eps1(nwdos),eps2(nwdos))

    ! generate energy grids
    brd=0.d0
    if (acont) brd=brdtd
    call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
    call genwgrid(nwdos,wdos,.false.,brd,w_cmplx=wr)
    wplot=dble(wr)

    if (nproc.gt.1) then
       write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, rank
    else
       write(filextp,'("_Q",i5.5,".OUT")') iq
    end if

    ! record length
    inquire(iolength=recl) mdf1(1)
    call getunit(unit1)

    ! neglect/include local field effects
    do m=1,n,max(n-1,1)

       ! loop over longitudinal components for optics
       do oct=1,nc

          optcomp(1,1)=oct
          optcomp(2,1)=oct
          optcompt(:)=optcomp(:,1)
          ! symmetrization matrix for dielectric function
          call gensymdf(oct,oct)

          ! string for xc-kernel type
          write(str,'(i2.2)') fxctype
          if (tq0) write(str,'(i2.2,"_OC",i2.2)') fxctype,11*oct
          str='_FXC'//trim(str)
          if (m.eq.1) str='_NLF'//trim(str)
          if (.not.aresdf) str='_NAR'//trim(str)
          if (acont) str='_AC'//trim(str)
          if (tetra) str='_TET'//trim(str)
          filnam2='IDF'//trim(str)

          ! read macroscopic dielectric function (original frequencies)
          open(unit1,file=trim(filnam2)//trim(filext),form='unformatted', &
               action='read',status='old',access='direct',recl=recl)
          do iw=1,nwdf
             read(unit1,rec=iw) mdf1(iw)
          end do
          close(unit1)

          ! analytic continuation
          if (acont) then
             call pade(nwdos,wr,nwdf,w,mdf1,mdf)
          else
             mdf(:)=mdf1(:)
          end if

          ! file names for optical functions
          fneps='EPSILON'//trim(str)
          fnloss='LOSS'//trim(str)
          fnsigma='SIGMA'//trim(str)
          fnsumrules='SUMRULES'//trim(str)
          if (tetra.and.(fxctype/=0)) then
             write(str,'(i2.2)') 0
             if (tq0) write(str,'(i2.2,"_OC",i2.2)') 0,11*oct
             str='_FXC'//trim(str)
             str='_NLF'//trim(str)
             if (.not.aresdf) str='_NAR'//trim(str)
             str='_TET'//trim(str)
             filnam='IDF'//trim(str)

             ! read macroscopic dielectric function (RPA)
             open(unit1,file=trim(filnam)//trim(filext),form='unformatted', &
                  action='read',status='old',access='direct',recl=recl)
             do iw=1,nwdf
                read(unit1,rec=iw) zt1
                mdfrpa(iw,oct,1)=dble(zt1)
                mdfrpa(iw,oct,2)=aimag(zt1)
             end do
             close(unit1)
             ! derivative of RPA dielectric function
             call fderiv(1,nwdos,wplot,mdfrpa(1,oct,1),mdfrpad(1,oct),cf)
             ! derivative of xc-kernel
             call fderiv(1,nwdos,wplot,fxc0(1,oct),fxc0d(1,oct),cf)

             write(str,'(i2.2)') fxctype
             if (tq0) write(str,'(i2.2,"_OC",i2.2)') fxctype,11*oct
             str='_FXC'//trim(str)
             str='_NLF'//trim(str)
             if (.not.aresdf) str='_NAR'//trim(str)
             str='_TET'//trim(str)
             fnexciton2='EXCITON'//trim(str)

             ! determination of excitons in combination with tetrahedron method
             ! and neglect local field effects for kernel but consider for
             ! the RPA solution
             if (m==max(n,1)) then
                call findexciton(oct,nwdos,dble(w))
                call writeexciton(iq,oct,wplot,mdf,trim(fnexciton2)//&
                     trim(filext))
             end if             
          end if

          ! generate optical functions
          call genloss(mdf,loss)
          call gensigma(dble(wr),mdf,optcompt,sigma)
          call gensumrls(dble(wr),mdf,sumrls)
          ! write optical functions to file
          call writeeps(iq,wplot,mdf,trim(fneps)//trim(filext))
          call writeloss(iq,wplot,loss,trim(fnloss)//trim(filext))
          call writesigma(iq,wplot,sigma,trim(fnsigma)//trim(filext))
          call writesumrls(iq,sumrls,trim(fnsumrules)//trim(filext))
       end do ! oct
    end do ! m
    ! deallocate
    if (allocated(mdfrpa)) deallocate(mdfrpa)
    deallocate(mdf,mdf1,w,wr,wplot,loss,sigma)
    deallocate(eps1,eps2,cf)

    ! restore file extension
    write(filext,'(".OUT")')

  end subroutine tdlinopt

end module m_tdlinopt
