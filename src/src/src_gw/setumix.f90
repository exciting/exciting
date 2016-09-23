!BOP
!
!!ROUTINE: setumix
!
!!INTERFACE:
!
subroutine setumix(ia,is)
!      
!!DESCRIPTION:
!
!Set up the radial part of the mixed basis functions ($\upsilon_{aNL}$)used for the matrix
!expansion of the non local operators (Polarization, Bare and Screened Coulomb
!potential, and the Self energy) for atom \verb"ias". The procedure is as follows:
!
!\begin{itemize}
!\item For each $L$ we take the product of radial functions $u_{lm}(r)u_{l'm'}(r)$ which fullfil the
!condition $|l-l'|\le L \le l+l'$.
!
!\item We calculate the overlap matrix of the products of radial functions:
!\begin{equation}
!\mathbb{O}_{(l,l');(l_1,l'_1)}=\int\limits_0^{R^a_{MT}} {u_{al}u_{al'}(\vec{r})u_{al_1}u_{al'_1}(\vec{r})
!r^2 dr}
!\end{equation}
!
!\item Diagonalize the matrix $\mathbb{O}_{(lm,l'm');(l_1m_1,l'_1m'_1)}$
!
!\item Discard the eigenvectors corresponding to eigenvalues with absolute value lower than a given
!tolerance (usually $10^{-5}$)
!
!\item The rest of the eigenvectors are normalized and stored for a grid of $\vec{r}$ that constitute the
!new basis $\left\{\upsilon_{NL}\right\}$
!
!\end{itemize}
!
!So defined the set of functions $\left\{\gamma
!_{aNLM}=\upsilon_{aNL}Y_{LM}\right\}$ constitute and orthonormal basis set,
!that is:
!
!\begin{equation}\label{mixborth}
!\int\limits_{V^a_{MT}}\gamma _{aNLM}(\vec{r})\gamma
!_{aN'L'M'}(\vec{r})d^3r=\delta_{N,N'}\delta_{L,L'}\delta_{M,M'}
!\end{equation}
!
!!USES:      
    use modinput
    use modmain
    use modmpi
    use modgw
    use reallocate

!!INPUT VARIABLES:
    implicit none
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: is
    
!!LOCAL VARIABLES:
    integer(4) :: ias
    integer(4) :: ir,i,j,k,l
    integer(4) :: nwf
    integer(4) :: nl_,nli,nor
    integer(4) :: lammax,lammin
    integer(4) :: i1,i2,info,j2
    integer(4) :: nl(0:2*maxlapw)
    real(8) :: fr(nrmtmax)
    real(8) :: gr(nrmtmax) 
    real(8) :: cf(3,nrmtmax)
   
    integer(4), allocatable :: ind(:)
    real(8),    allocatable :: uol(:), work(:)
    real(8),    allocatable :: uml(:,:)

!!EXTERNAL ROUTINES: 
    external dsyev
      
!!REVISION HISTORY:
!
! Created May. 19th. 2004 by RGA
! Last Modified May 25th 2004 by RGA
! Revisited 28.04.2011 by DIN
! Modified Dec 2013 by DIN
!
!EOP
!BOC

    ias = idxas(ia,is)
        
    ! Calculate the number of radial product functions for each L 
    ! and the total number of orbitals
    nl(:) = 0
    do i= 1, nup
      lammax = eles(i,1)+eles(i,2)
      lammin = abs(eles(i,1)-eles(i,2))
      do l = 0, 2*input%groundstate%lmaxapw
        if ((lammin<=l).and.(l<=lammax)) then
          nl(l) = nl(l)+1
        end if 
      end do ! l
    end do
          
    !-----------------------------------------------------------------------
    ! Diagonalize the overlap matrix of the product functions for each L-block
    !-----------------------------------------------------------------------
    j2 = 0
    nor = 0
    do l = 0, 2*input%groundstate%lmaxapw
        
      nl_ = nl(l)
        
      if (nl_>1) then

        ! allocate the L-block overlap matrix (uml), the eigenvalues
        ! vector (uol), the working space and the indexes
        allocate(uml(nl_,nl_))
        uml(:,:) = 0
        allocate(uol(nl_))
        uol(:) = 0
        allocate(ind(nl_))
        ind(:) = 0
              
        ! generate the L-block overlap matrix
        j = 0
        do i1 = 1, nup
          lammax = eles(i1,1)+eles(i1,2)
          lammin = abs(eles(i1,1)-eles(i1,2))
          if ((lammin<=l).and.(l<=lammax)) then
            j = j+1
            ind(j) = i1
            k = j
            uml(j,j) = umat(i1,i1)
            do i2 = i1+1, nup
              lammax = eles(i2,1)+eles(i2,2)
              lammin = abs(eles(i2,1)-eles(i2,2))
              if ((lammin<=l).and.(l<=lammax)) then
                k = k+1
                uml(j,k) = umat(i1,i2)
                uml(k,j) = uml(j,k)
              end if
            end do ! i2 
          end if
        end do ! i1
              
        if (input%gw%debug) then
          write(fdebug,*) "###uml for L=", l
          do k = 1, nl_
            write(fdebug,'(100f16.8)') uml(1:nl_,k)
          end do
        end if
             
        ! diagonalize the L-block overlap matrix
        allocate(work(3*nl_-1))
        call dsyev('V','U',nl_,uml,nl_,uol,work,3*nl_-1,info)
        if (info.ne.0) then
          write(*,*) 'setumix: error in calling dsyev '
          write(*,*) 'l =', l,' nl =', nl_, 'info = ',info
          stop              
        end if

#ifdef MPI
        ! (chm): the result uml can vary because math libs can be nondeterministic
        ! this broadcast just ensures that each process starts with the exact same data
        call MPI_BCAST(uml,nl_*nl_,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
            
        !---------------------------------------
        ! Get rid of linear dependent functions
        !---------------------------------------
        i2 = 0
        do i1 = 1, nl_
          if (uol(i1)>input%gw%MixBasis%epsmb) i2 = i2+1
        end do 
        nli = i2
        
        if (input%gw%debug) write(fdebug,101) l, nl_, nli

        ! transform the umix 
        do i1 = 1, nl_
          if (uol(i1)>input%gw%MixBasis%epsmb) then
            j2 = j2+1
            do i2 = 1, nl_
              do ir = 1, nrmt(is)
                umix(ir,j2,ias) = umix(ir,j2,ias)+ &
                &                 uml(i2,i1)*uprod(ir,ind(i2))
              end do ! ir
            end do ! i2 
          end if
        end do ! i1

        !--------------------------------------------
        ! normalize the radial mixed basis functions
        !--------------------------------------------
        do i1 = nor+1, nor+nli
          bigl(i1,ias) = l
          do ir = 1, nrmt(is)
            fr(ir) = umix(ir,i1,ias)*umix(ir,i1,ias)
          end do ! ir
          call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
          umix(1:nrmt(is),i1,ias) = umix(1:nrmt(is),i1,ias)/sqrt(gr(nrmt(is)))
        end do ! i1
            
        if (input%gw%debug) then
          do ir = 1, nrmt(is)
            write(fdebug,*) spr(ir,is), (umix(ir,i1,ias)/spr(ir,is),i1=nor+1,nor+nli)
          end do
        end if
              
        nor = nor+nli
        nl(l) = nli

        ! deallocate the local arrays
        deallocate(uml)
        deallocate(uol)
        deallocate(work)
        deallocate(ind)

      !--------------------------------------------------
      ! in case the L-block is just one basis functions
      !--------------------------------------------------
      else if (nl_==1) then
              
        if (input%gw%debug) write(fdebug,101) l, nl_
              
        do i1 = 1, nup
          lammax = eles(i1,1)+eles(i1,2)
          lammin = abs(eles(i1,1)-eles(i1,2))
          if((lammin<=l).and.(l<=lammax))then
            j2 = j2+1
            bigl(j2,ias) = l
            do ir = 1, nrmt(is)
              fr(ir) = uprod(ir,i1)*uprod(ir,i1)
            end do ! ir
            call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
            umix(1:nrmt(is),j2,ias) = uprod(1:nrmt(is),i1)/sqrt(gr(nrmt(is)))
          end if
        end do ! i1
        nor = nor+1
            
        if (input%gw%debug) then
          write(fdebug,*)
          do ir = 1, nrmt(is)
            write(fdebug,10) spr(ir,is), umix(ir,j2,ias)/spr(ir,is)
          end do
        end if
              
      end if ! nl>1
    enddo ! l

    ! Store in the global arrays
    nmix(ias) = nor
    mbl(ias) = maxval(bigl(:,ias))
    
    ! Deallocate not needed global arrays
    deallocate(uprod)
    deallocate(umat)
    deallocate(eles)
    
    ! print debugging info
    if (input%gw%debug) then
      write(fdebug,102) nmix(ias), mbl(ias)
      write(fdebug,103)
      nwf = 0
      do i = 1, nor
        write(fdebug,104) i, bigl(i,ias), 2*bigl(i,ias)+1
        nwf = nwf+2*bigl(i,ias)+1
      end do
      write(fdebug,105) nwf
    end if    
      
    10  format(41d15.7)
    100 format(/,10x,'Mixed basis functions for atom',1x,a10,/)
    101 format(10x,'L =',i2,' Nr. of products:',i4,' Nr. of basis functions',i4)  
    102 format(59x,'----',/,10x,'Total number of radial functions',17x,   &
    &          i4,4x,'Maximum L',i4)
    103 format(26x,'N',3x,'L',2x,'deg.')
    104 format(23x,3i4)
    105 format(31x,'----',/,'Total number of basis functions',i4)
    
    return
end subroutine
!EOC
