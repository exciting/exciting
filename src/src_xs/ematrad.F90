
module m_ematrad
  implicit none
contains

  subroutine ematrad(iq)
    use modmain
    use modtddft
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    integer is,ia,ias,nr,ir,igq
    integer l1,l2,l3
    integer ilo,ilo1,ilo2,io,io1,io2
    real(8) t1
    integer :: lmax1, lmax2, lmax3
    integer :: u11,u22,u33
    ! automatic arrays
    real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
    ! allocatable arrays
    real(8), allocatable :: jl(:,:), jhelp(:)

    lmax1 = lmaxapwtd
    lmax2 = lmaxemat
    ! lmax1 and lmax3 should be the same!
    lmax3 = lmax1

    ! allocate arrays for radial integrals and Bessel functions
    if (allocated(riaa)) deallocate(riaa)
    if (allocated(riloa)) deallocate(riloa)
    if (allocated(rilolo)) deallocate(rilolo)
    allocate(riaa(0:lmax1,apwordmax,0:lmax3,apwordmax,0:lmax2,natmtot,ngq(iq)))
    allocate(riloa(nlomax,0:lmax3,apwordmax,0:lmax2,natmtot,ngq(iq)))
    allocate(rilolo(nlomax,nlomax,0:lmax2,natmtot,ngq(iq)))
    ! allocate temporary arrays
    allocate(jl(0:lmax2,nrmtmax))
    allocate(jhelp(0:lmax2))
    jl(:,:) = 0.d0
    jhelp(:) = 0.d0
    ! zero arrays for radial integrals
    riaa(:,:,:,:,:,:,:) = 0.d0
    riloa(:,:,:,:,:,:) = 0.d0
    rilolo(:,:,:,:,:) = 0.d0

    if (dbglev.gt.1) then
       ! APW-APW
       call getunit(u11)
       open(unit=u11,file='IRADaa'//filext,form='formatted',action='write', &
            status='replace')
       write(u11,'(a)') 'igq,ias,l1,io1,l3,io2,l2   iraa'
       write(u11,'(a)') '-----------------------------------------------------'
       ! lo-APW
       call getunit(u22)
       open(unit=u22,file='IRADalo'//filext,form='formatted',action='write', &
            status='replace')
       write(u22,'(a)') 'igq,ias,ilo,l1,l3,io,l2,   iralo'
       write(u22,'(a)') '-----------------------------------------------------'
       ! lo-lo
       call getunit(u33)
       open(unit=u33,file='IRADlolo'//filext,form='formatted',action='write', &
            status='replace')
       write(u33,'(a)') 'igq,ias,ilo1,l1,ilo2,l3,l2,   irlolo'
       write(u33,'(a)') '-----------------------------------------------------'
    end if

    ! begin loop over G+q vectors
    do igq=1,ngq(iq)
       ! begin loop over species
       do is=1,nspecies
          nr=nrmt(is)
          do ir=1,nr
             ! calculate r^2
             r2(ir)=spr(ir,is)**2
             ! calculate spherical Bessel functions of first kind j_l(|G+q|r_a)
             call sbessel(lmax2,gqc(igq,iq)*spr(ir,is),jhelp)
             jl(:,ir) = jhelp(:)
          end do
          ! begin loop over atoms
          do ia=1,natoms(is)
             ias=idxas(ia,is)
             !----------------!
             !     APW-APW    !
             !----------------!
             do l1=0,lmax1
                do io1=1,apword(l1,is)
                   do l3=0,lmax3
                      do io2=1,apword(l3,is)
                         do l2=0,lmax2
                            do ir=1,nr
                               t1=apwfr(ir,1,io1,l1,ias)* &
                                    apwfr(ir,1,io2,l3,ias)*r2(ir)
                               fr(ir)=t1*jl(l2,ir)
                            end do
                            call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                            riaa(l1,io1,l3,io2,l2,ias,igq)=gr(nr)
                            !irad(io1,l2,l3,io2)=gr(nr)
                            if (dbglev.gt.1) then
                               write(u11,'(7i5,g18.10)') &
                                    igq,ias,l1,io1,l3,io2,l2, &
                                    gr(nr)
                            end if
                         end do
                      end do ! io2
                   end do ! l3
                end do ! io1
             end do ! l1
             !----------------------------!
             !     local-orbital-APW      !
             !----------------------------!
             do ilo=1,nlorb(is)
                l1=lorbl(ilo,is)
                do l3=0,lmax3
                   do io=1,apword(l3,is)
                      do l2=0,lmax2
                         do ir=1,nr
                            t1=lofr(ir,1,ilo,ias)* &
                                 apwfr(ir,1,io,l3,ias)*r2(ir)
                            fr(ir)=t1*jl(l2,ir)
                         end do
                         call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                         riloa(ilo,l3,io,l2,ias,igq)=gr(nr)
                         !irad(1,l2,l3,io)=gr(nr)
                         if (dbglev.gt.1) then
                            write(u22,'(7i5,g18.10)') &
                                 igq,ias,ilo,l1,l3,io,l2, &
                                 gr(nr)
                         end if
                      end do ! l2
                   end do ! io
                end do ! l3
             end do ! ilo
             !------------------------------------!
             !     local-orbital-local-orbital    !
             !------------------------------------!
             do ilo1=1,nlorb(is)
                l1=lorbl(ilo1,is)
                do ilo2=1,nlorb(is)
                   l3=lorbl(ilo2,is)
                   do l2=0,lmax2
                      do ir=1,nr
                         t1=lofr(ir,1,ilo1,ias)* &
                              lofr(ir,1,ilo2,ias)*r2(ir)
                         fr(ir)=t1*jl(l2,ir)
                      end do
                      call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                      rilolo(ilo1,ilo2,l2,ias,igq)=gr(nr)
                      !irad(1,l2,l3,ilo2)=gr(nr)
                      if (dbglev.gt.1) then
                         write(u33,'(7i5,g18.10)') &
                              igq,ias,ilo1,l1,ilo2,l3,l2, &
                              gr(nr)
                      end if
                   end do ! l2
                end do ! ilo2
             end do ! ilo1
             ! end loops over atoms and species
          end do
       end do
       ! end loop over G+q vectors
    end do

    ! deallocate
    deallocate(jl,jhelp)
    if (dbglev.gt.1) then
       ! close files
       close(u11)
       close(u22)
       close(u33)
    end if

  end subroutine ematrad

end module m_ematrad
