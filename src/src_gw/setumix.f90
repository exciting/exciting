!BOP
!
! !ROUTINE: setumix
!
! !INTERFACE:
      subroutine setumix()
      
! !DESCRIPTION:
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
!

! !USES:      

      use modinput
      use modmain
      use modgw
      use modmpi
      use reallocate
 
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: ir
      integer(4) :: is
      integer(4) :: jprod
      integer(4) :: nuor
      integer(4) :: i,j,k,l
      integer(4) :: nli,nor
      integer(4) :: lammax, lammin, nwf
      integer(4) :: i1,i2,info,j2

      integer(4) :: nl(0:2*maxlapw)
      real(8) :: fr(nrmtmax)
      real(8) :: gr(nrmtmax) 
      real(8) :: cf(3,nrmtmax)
      
      integer(4), allocatable :: ind(:)
      real(8), allocatable :: uol(:),work(:)
      real(8), allocatable :: uml(:,:)

!
! !EXTERNAL ROUTINES: 

      external dsyev
      external fderiv
      external ylm
      external setuprod
      
!
! !REVISION HISTORY:
!
! Created May. 19th. 2004 by RGA
! Last Modified May 25th 2004 by RGA
! Revisited 28.04.2011 by DIN
!
!EOP
!BOC
!
      if (allocated(nmix)) deallocate(nmix)
      allocate(nmix(natmtot))
      if (associated(umix)) deallocate(umix)
      allocate(umix(natmtot,maxnup,nrmtmax))
      if (associated(bigl)) deallocate(bigl)
      allocate(bigl(natmtot,maxnup))
      if (allocated(mbl)) deallocate(mbl)
      allocate(mbl(natmtot))
      maxnmix=0
      maxbigl=0

      do is=1,nspecies
        if(debug)write(701,100)spname(is)
        do ia=1,natoms(is)
          ias=idxas(ia,is)
!
!     Calculate the number of radial product functions for each L and the 
!     total number of orbitals
!
          nmix(ias)=0
          nl(:)=0
          do jprod=1,nup(ias)
            lammax=eles(ias,jprod,1)+eles(ias,jprod,2)
            lammin=abs(eles(ias,jprod,1)-eles(ias,jprod,2))
            do l=0,2*input%groundstate%lmaxapw
              if((l.le.lammax).and.(l.ge.lammin))then
                nl(l)=nl(l)+1
              endif 
            enddo ! l
          enddo ! jprod

          umix(ias,:,:)=0.0d0
          bigl(ias,:)=0
!
!     diagonalize the overlap matrix of the product functions for each  L-block
!
          j2=0
          nor=0
          do l=0,2*input%groundstate%lmaxapw
            if (nl(l).gt.1) then
!
!             allocate the L-block overlap matrix (uml), the eigenvalues
!             vector (uol), the working space and the indexes
!
              allocate(uml(nl(l),nl(l)))
              uml=0
              allocate(uol(nl(l)))
              uol=0
              allocate(work(3*nl(l)-1))
              allocate(ind(nl(l)))
!
!             generate the L-block overlap matrix
!
				 
              j=0
              do i1=1,nup(ias)
                lammax=eles(ias,i1,1)+eles(ias,i1,2)
                lammin=abs(eles(ias,i1,1)-eles(ias,i1,2))
                if((l.le.lammax).and.(l.ge.lammin))then
                  j=j+1
                  ind(j)=i1
                  k=j
                  uml(j,j)=umat(ias,i1,i1)
                  do i2=i1+1,nup(ias)
                    lammax=eles(ias,i2,1)+eles(ias,i2,2)
                    lammin=abs(eles(ias,i2,1)-eles(ias,i2,2))
                    if((l.le.lammax).and.(l.ge.lammin))then
                      k=k+1
                      uml(j,k)=umat(ias,i1,i2)
                      uml(k,j)=uml(j,k)
                    endif
                  enddo ! i2 
                endif
              enddo ! i1
              
              if(debug)then
                write(701,*) "###uml for L=", l
                do k=1,nl(l)
                  write(701,'(100f16.8)') uml(1:nl(l),k)
                enddo
              end if
              
!             diagonalize the L-block overlap matrix
 
			 
              call dsyev('V','U',nl(l),uml,nl(l),uol,work,3*nl(l)-1,info)
             
              if(info.ne.0)then
                write(*,*) 'setumix: error in calling dsyev '
                write(*,*) 'l =',l,' nl =',nl(l),'info = ',info
                stop
              endif
#ifdef MPI
			! the result uml can vary because math libs can be nondeterministic
			! this broadcast just ensures that each process starts with the exact same data
             call mpi_bcast(uml,nl(l)*nl(l),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
              
              i2=0
              do i1=1,nl(l)
                if (uol(i1).gt.input%gw%MixBasis%epsmb) i2=i2+1
              enddo 
              nli=i2
              
              if(debug)write(701,101) l,nl(l),nli
!
!             transform the umix 
!
              do i1=1,nl(l)
                if(uol(i1).gt.input%gw%MixBasis%epsmb)then
                  j2=j2+1
                  do i2=1,nl(l)
                    do ir=1,nrmt(is)
                      umix(ias,j2,ir)=umix(ias,j2,ir)+uml(i2,i1)*uprod(ias,ind(i2),ir)
                    enddo ! ir
                  enddo ! i2 
                endif
              enddo ! i1
           
              
!
!             normalize the radial mixed wave functions
!
              do i1=nor+1,nor+nli
                bigl(ias,i1)=l
                do ir=1,nrmt(is)
                   fr(ir)=umix(ias,i1,ir)*umix(ias,i1,ir)
                enddo ! ir
                call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                umix(ias,i1,1:nrmt(is))=umix(ias,i1,1:nrmt(is))/sqrt(gr(nrmt(is)))
              enddo ! i1
            
!             write the mixed wave functions to disk          
              if(debug)then
                do ir=1,nrmt(is)
                  write(701,*)spr(ir,is),(umix(ias,i1,ir)/spr(ir,is),i1=nor+1,nor+nli)
                enddo
              end if
              
              nor=nor+nli
              nl(l)=nli
!
!             deallocate the temporary arrays
!
              deallocate(uml)
              deallocate(uol)
              deallocate(work)
              deallocate(ind)
!
!           in case the L-block is just one wavefunctions
!
            else if (nl(l).eq.1) then
              
              if(debug)write(701,101) l,nl(l)
              
              do i1=1,nup(ias)
                lammax=eles(ias,i1,1)+eles(ias,i1,2)
                lammin=abs(eles(ias,i1,1)-eles(ias,i1,2))
                if((l.le.lammax).and.(l.ge.lammin))then
                  j2=j2+1
                  bigl(ias,j2)=l
                  do ir=1,nrmt(is)
                    fr(ir)=uprod(ias,i1,ir)*uprod(ias,i1,ir)
                  enddo ! ir
                  call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                  umix(ias,j2,:)=uprod(ias,i1,:)/sqrt(gr(nrmt(is)))
                endif
              enddo ! i1
              nor=nor+1
             
!             write the wave function to disk
              if(debug)then
                write(701,*)  
                do ir=1,nrmt(is)
                  write(701,11)spr(ir,is),umix(ias,j2,ir)/spr(ir,is)
                enddo
              end if
              
            endif
          enddo ! l

          nmix(ias)=nor
          if(nmix(ias).gt.maxnmix) maxnmix=nmix(ias)
          
          mbl(ias)=maxval(bigl(ias,:))
          if(mbl(ias).gt.maxbigl) maxbigl=mbl(ias)
          
          if(debug)then
            write(701,102) nmix(ias), mbl(ias)
            write(701,103)
            nwf=0
            do i=1,nor
              write(701,104)i,bigl(ias,i),2*bigl(ias,i)+1
              nwf=nwf+2*bigl(ias,i)+1
            enddo
            write(701,105) nwf
          end if

        enddo ! ia
      enddo ! is
      
      call doreallocate(umix,natmtot,maxnmix,nrmtmax)
      call doreallocate(bigl,natmtot,maxnmix)
      
      deallocate(uprod)
      deallocate(umat)
      deallocate(eles)
      
   11 format(41d15.7)
  100 format(/,10x,'Mixed basis functions for atom',1x,a10,/)
  101 format(10x,'L =',i2,' Nr. of products:',i4,' Nr. of basis functions',i4)  
  102 format(59x,'----',/,10x,'Total number of radial functions',17x,   &
     &       i4,4x,'Maximum L',i4)
  103 format(26x,'N',3x,'L',2x,'deg.')
  104 format(23x,3i4)
  105 format(31x,'----',/,'Total number of basis functions',i4)
      return

      end subroutine setumix
!EOC
