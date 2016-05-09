
!--------------------------------------------!
!     Brillouin zone integration             !
!--------------------------------------------!

module mod_bzintw
    
!-------------------------------------------------------------------------------
    type k_bzintw
    
        !------------------
        ! Core electrons
        !------------------
        ! k-dependent integration weight used for the summation like 
        ! $\sum_k^{\mathrm{BZ}} \sum_n^{\mathrm{occ}} F_{nk} $
        real(8), allocatable :: ciw(:,:)

        !------------------
        ! Valence electrons
        !------------------
        ! k-dependent integration weight used for the summation like 
        ! $\sum_k^{\mathrm{BZ}} \sum_n^{\mathrm{occ}} F_{nk} $
        real(8), allocatable :: kiw(:,:)
    
        ! k-point weight of a certain band for a normal surface integration   
        real(8), allocatable :: kwfer(:,:)

    end type k_bzintw
    
!-------------------------------------------------------------------------------
    type kq_bzintw
    
        !------------------
        ! Core electrons
        !------------------
        ! k/q-dependent integration weights used for the summation like
        ! $\sum_k^{\mathrm{BZ}} \sum_n^{\mathrm{occ}} \sum_m^{\mathrm{unocc}} F_{nk,mk-q} $
        complex(8), allocatable :: unw(:,:,:,:,:)   

        !------------------
        ! Valence electrons
        !------------------
        ! k/q-dependent integration weights used for the summation like
        ! $\sum_k^{\mathrm{BZ}} \sum_n^{\mathrm{occ}} \sum_m^{\mathrm{unocc}} F_{nk,mk-q} $
        complex(8), allocatable :: kcw(:,:,:,:) 
    
    end type kq_bzintw
    
contains

!-------------------------------------------------------------------------------    
    subroutine generate_k_bzintw(self,kset,lcore)
        use modmain, only: nspecies, natoms, idxas, &
        &                  spnstmax, spnst, spcore, natmtot, &
        &                  evalcr, nstsv, evalsv, efermi
        use mod_kpointset
        implicit none
        type(k_bzintw), intent(out) :: self
        type(k_set),    intent(in)  :: kset
        logical,        intent(in)  :: lcore
        ! local variables        
        integer :: ia, is, ias
        integer :: ic, ist, ik, ikp
        real(8), allocatable :: cwpar(:,:)
        real(8), allocatable :: bandpar(:,:)
        
        !-------------------
        ! core
        !-------------------
        if (lcore) then
          if (allocated(self%ciw)) deallocate(self%ciw)
          allocate(self%ciw(natmtot,spnstmax))
          self%ciw(:,:) = 0.0d0
          allocate(bandpar(1,kset%nkpt))
          allocate(cwpar(1,kset%nkpt))
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              ic = 0
              do ist = 1, spnst(is)
                if (spcore(ist,is)) then
                  ic = ic+1
                  bandpar(1,:) = evalcr(ic,ias)
                  call tetiw(     &
                  &  kset%nkpt,   &
                  &  kset%ntet,   &
                  &  1, bandpar,  &
                  &  kset%tnodes, &
                  &  kset%wtet,   &
                  &  kset%tvol,   &
                  &  efermi,      &
                  &  cwpar )
                  self%ciw(ias,ic) = cwpar(1,1)
                end if
              enddo ! ist
            enddo ! ia
          enddo ! is
          deallocate(bandpar)
          deallocate(cwpar)
        end if ! core

        !--------------------
        ! valence
        !--------------------
        if (allocated(self%kiw)) deallocate(self%kiw)
        allocate(self%kiw(nstsv,kset%nkpt))
        self%kiw = 0.0d0
        if (allocated(self%kwfer)) deallocate(self%kwfer)        
        allocate(self%kwfer(nstsv,kset%nkpt))
        self%kwfer = 0.0d0
        allocate(bandpar(nstsv,kset%nkpt))
        do ik = 1, kset%nkpt
          ikp = kset%ik2ikp(ik)
          bandpar(1:nstsv,ik) = evalsv(1:nstsv,ikp)
        end do
        call tetiw( &
        &  kset%nkpt, &
        &  kset%ntet, &
        &  nstsv, bandpar, &
        &  kset%tnodes, &
        &  kset%wtet, &
        &  kset%tvol, &
        &  efermi, &
        &  self%kiw )
        call tetiwsurf( &
        &  kset%nkpt, &
        &  kset%ntet, &
        &  nstsv, bandpar, &
        &  kset%tnodes, &
        &  kset%wtet, &
        &  kset%tvol, &
        &  efermi, &
        &  self%kwfer )
        deallocate(bandpar)
        
        return
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine delete_k_bzintw(self)
        type(k_bzintw), intent(inout) :: self
        if (allocated(self%ciw)) deallocate(self%ciw)
        if (allocated(self%kiw)) deallocate(self%kiw)
        if (allocated(self%kwfer)) deallocate(self%kwfer)
    end subroutine

!-------------------------------------------------------------------------------

    subroutine print_k_bzintw(self,funit)
        use modmain, only: nspecies, natoms, idxas, &
        &                  spnst, spcore, nstsv 
        implicit none
        type(k_bzintw), intent(in) :: self
        integer,        intent(in) :: funit
        ! local varible
        integer :: ia, is, ias, nk, ik, ist, ic
        
        call boxmsg(funit,'-','BZ integration weights')
        if (allocated(self%ciw)) then
          call linmsg(funit,'-','Core states: < atom    state     ciw >')
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              ic = 0
              do ist = 1, spnst(is)
                if (spcore(ist,is)) then
                  ic = ic+1
                  write(funit,'(2i6,f12.6)') ias, ist, self%ciw(ias,ic)
                end if
              end do ! ist
            end do ! ia
          end do ! is
        endif
        call linmsg(funit,'-','Valence states: < ik    state     kiw >')
        nk = size(self%kiw,2)
        do ik = 1, nk
          do ist = 1, nstsv
            write(funit,'(2i6,f12.6)') ik, ist, self%kiw(ist,ik)
          end do ! ist
        end do ! ik
        return
    end subroutine

!-------------------------------------------------------------------------------
    subroutine generate_kq_bzintw(self,iq,kset,kqset,freq,nomax,numin,lcore)
        use modmain, only: nspecies, natoms, idxas, spnstmax, & 
        &                  spnst, spcore, natmtot, &
        &                  evalcr, nstsv, evalsv, efermi
        use mod_kpointset
        use mod_frequency
        implicit none
        type(kq_bzintw), intent(out) :: self
        integer,         intent(in)  :: iq
        type(k_set),     intent(in)  :: kset
        type(kq_set),    intent(in)  :: kqset
        type(frequency), intent(in)  :: freq
        integer,         intent(in)  :: nomax
        integer,         intent(in)  :: numin
        logical,         intent(in)  :: lcore
        ! local variables        
        integer :: ist, jst, ic, iom
        integer :: fflg
        integer :: ik, ikp, ia, is, ias, sgw
        real(8) :: emaxb, edif, edsq, omsq
        real(8), allocatable :: cwpar(:,:,:)
        real(8), allocatable :: cwparsurf(:,:,:)
        real(8), allocatable :: bandpar(:,:)

        ! Initialization        
        select case (freq%fconv)
          case('nofreq')
            fflg = 1
          case('refreq')
            fflg = 2
          case('imfreq')
            fflg = 3
        end select
        sgw = 5-2*fflg
        
        write(*,*) freq%nomeg
        write(*,*) kqset%nkpt

        !---------------
        ! Core-Valence     
        !---------------
        if (lcore) then 
          if (allocated(self%unw)) deallocate(self%unw)
          allocate(self%unw(natmtot,spnstmax,nstsv,freq%nomeg,kqset%nkpt))
          self%unw(:,:,:,:,:) = 0.d0
          allocate(bandpar(2,kqset%nkpt))
          allocate(cwpar(2,2,kqset%nkpt))
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              ic = 0
              do ist = 1, spnst(is)
                if (spcore(ist,is)) then
                  ic = ic+1
                  bandpar(1,:) = evalcr(ic,ias)
                  do jst = numin, nstsv
                    do ik = 1, kqset%nkpt
                      ikp = kset%ik2ikp(ik)
                      bandpar(2,ik) = evalsv(jst,ikp)
                    end do
                    ! continue only if the band is at least partially unoccupied
                    emaxb = maxval(bandpar(2,:))
                    if (emaxb > efermi) then
                      do iom = 1, freq%nomeg
                        omsq = sgw*freq%freqs(iom)*freq%freqs(iom)
                        ! libbzint routine
                        ! old version
                        !call tetcw( &
                        !&  kqset%nkpt, &
                        !&  kqset%ntet, &
                        !&  2, &
                        !&  kqset%wtet, &
                        !&  bandpar, &
                        !&  kqset%tnodes, &
                        !&  kqset%linkq(:,iq), &
                        !&  kqset%tvol, &
                        !&  efermi,&
                        !&  freq%freqs(iom), &
                        !&  1, cwpar )
                        ! new version
                        call tetcw( &
                        &  kqset%nkpt, &
                        &  kqset%ntet, &
                        &  2, &
                        &  kqset%wtet, &
                        &  bandpar, &
                        &  kqset%tnodes, &
                        &  kqset%linkq(:,iq), &
                        &  kqset%kqid(:,iq), &
                        &  kqset%tvol, &
                        &  efermi,&
                        &  freq%freqs(iom), &
                        &  1, cwpar )
     
                        do ik = 1, kqset%nkpt
                          ikp = kset%ik2ikp(ik)
                          edif = evalsv(jst,ikp)-evalcr(ic,ias)
                          edsq = edif*edif
                          self%unw(ias,ic,jst,iom,ik) = &
                          &  2.0d0*cwpar(1,2,ik)*edif/(omsq-edsq)
                        enddo
                        ! real frequencies                      
                        if (fflg==2) then
                          ! libbzint routine                       
                          call tetcw( &
                          &  kqset%nkpt, &
                          &  kqset%ntet, &
                          &  2, &
                          &  kqset%wtet, &
                          &  bandpar, &
                          &  kqset%tnodes, &
                          &  kqset%linkq(:,iq), &
                          &  kqset%kqid(:,iq), &
                          &  kqset%tvol, &
                          &  efermi,&
                          &  freq%freqs(iom), &
                          &  4, cwpar )
                       
                          do ik = 1, kqset%nkpt
                            self%unw(ic,ias,jst,iom,ik) = cmplx( &
                            &  real(self%unw(ic,ias,jst,iom,ik)), &
                            &  cwpar(1,2,ik))
                          enddo
                        endif
                      enddo ! iom
                    end if ! efermi
                  enddo ! jb
                end if
              enddo ! ist  
            enddo ! ia  
          enddo ! is
          deallocate(bandpar)
          deallocate(cwpar)
        end if ! core  

        !-------------------
        ! Valence-Valence     
        !-------------------
        if (allocated(self%kcw)) deallocate(self%kcw)
        allocate(self%kcw(nstsv,nstsv,1:freq%nomeg,kqset%nkpt))
        self%kcw(:,:,:,:) = 0.d0
        allocate(bandpar(nstsv,kqset%nkpt))
        do ik = 1, kqset%nkpt
          ikp = kset%ik2ikp(ik)
          bandpar(:,ik) = evalsv(:,ikp)
        end do
        allocate(cwpar(nstsv,nstsv,kqset%nkpt))
        if (fflg==2) then
          allocate(cwparsurf(nstsv,nstsv,kqset%nkpt))
        endif
        do iom = 1, freq%nomeg
          ! libbzint routine
          call tetcw( &
          &  kqset%nkpt, &
          &  kqset%ntet, &
          &  nstsv, &
          &  kqset%wtet, &
          &  bandpar, &
          &  kqset%tnodes, &
          &  kqset%linkq(:,iq), &
          &  kqset%kqid(:,iq), &
          &  kqset%tvol, &
          &  efermi,&
          &  freq%freqs(iom), &
          &  fflg, cwpar )
          ! real frequencies         
          if (fflg==2) then
            ! libbzint routine          
            call tetcw( &
            &  kqset%nkpt, &
            &  kqset%ntet, &
            &  nstsv, &
            &  kqset%wtet, &
            &  bandpar, &
            &  kqset%tnodes, &
            &  kqset%linkq(:,iq), &
            &  kqset%kqid(:,iq), &
            &  kqset%tvol, &
            &  efermi,&
            &  freq%freqs(iom), &
            &  4, cwparsurf )
          endif
          do ik = 1, kqset%nkpt
            do ist = 1, nstsv
              if (bandpar(ist,kqset%kqid(ik,iq)) > 900.0) &
              &  cwpar(:,ist,ik) = 0.0d0
              ! real frequencies             
              if (fflg==2) then
                if (bandpar(ist,kqset%kqid(ik,iq)) > 900.0) &
                &  cwparsurf(:,ist,ik) = 0.0d0
              endif
            enddo
          enddo
          ! real frequencies          
          if (fflg==2) then 
            self%kcw(1:nomax,numin:nstsv,iom,1:kqset%nkpt) = &
            &  cmplx(cwpar(1:nomax,numin:nstsv,1:kqset%nkpt), &
            &  cwparsurf(1:nomax,numin:nstsv,1:kqset%nkpt)) 
          else 
            self%kcw(1:nomax,numin:nstsv,iom,1:kqset%nkpt) = &
            &  cmplx(cwpar(1:nomax,numin:nstsv,1:kqset%nkpt),0.0)    
          endif 
        enddo ! iom
        deallocate(bandpar)
        deallocate(cwpar)
        if (fflg==2) deallocate(cwparsurf)
        return
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine delete_kq_bzintw(self)
        type(kq_bzintw), intent(inout) :: self
        if (allocated(self%unw)) deallocate(self%unw)
        if (allocated(self%kcw)) deallocate(self%kcw)
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine print_kq_bzintw(self,iq,funit)
        use modmain, only: nspecies, natoms, idxas, spnst, spcore, nstsv
        implicit none
        type(kq_bzintw), intent(in) :: self
        integer,         intent(in) :: iq
        integer,         intent(in) :: funit
        ! local varible
        integer :: ia, is, ias, no, nk, ik, ist, jst, ic, iom
        
        nk = size(self%kcw,4)
        no = size(self%kcw,3)
        call boxmsg(funit,'-','q-dependent BZ integration weights')
        write(funit,*) 'q-point < iq > : ', iq
        if (allocated(self%unw)) then
          call linmsg(funit,'-', &
          &  'Core-Valence: < iom    atom    core state    ik    val state    unw >')
          do ik = 1, nk
            do iom = 1, no
              do is = 1, nspecies
                do ia = 1, natoms(is)
                  ias = idxas(ia,is)
                  ic = 0
                  do ist = 1, spnst(is)
                    if (spcore(ist,is)) then
                      ic = ic+1
                      do jst = 1, nstsv
                        write(funit,'(5i6,2f12.6)') &
                        &  ik, iom, ias, ic, jst, self%unw(ias,ic,jst,iom,ik)
                      end do ! jst
                    end if
                  end do ! ist
                end do ! ia
              end do ! is
            end do ! iom
          end do ! ik
        endif
        call linmsg(funit,'-', &
        &  'Valence-Valence: < iom    val state    val state    ik    kiw >')
        do iom = 1, no
          do ist = 1, nstsv
            do jst = 1, nstsv
              do ik = 1, nk
                write(funit,'(4i6,2f12.6)') &
                &  iom, ist, jst, ik, self%kcw(ist,jst,iom,ik)
              end do ! ik
            end do ! jst
          end do ! ist
        end do ! iom
        return
    end subroutine

end module


