subroutine checkmbrot(iq)
    
    use modmain
    use modgw
    
    implicit none
    integer(4) :: iq, iqp
    integer(8) :: Recl
    integer(4) :: ik, ikp
    integer(4) :: isym, i
    integer(4) :: ist, jst, ijst, icg, dimtk
    
    complex(8), allocatable :: mmatk(:,:), mmatkp(:,:), rmmat(:,:), tmat(:,:)
    complex(8), allocatable :: m2k(:,:), m2r(:,:), m2tot(:,:)
    
    character(4) :: word
    
    ! irreducible q-point index
    iqp=indkpq(iq,1)

!--------------------------------------------
!   Check valence-valence M^i_{nm} rotation
!--------------------------------------------

    inquire(IoLength=Recl) minmmat
    open(41,file='minmmat.io',action='READ',form='UNFORMATTED', &
   &     access='DIRECT',recl=Recl)

    if(allocated(minmmat))deallocate(minmmat)
    allocate(minmmat(matsiz,nstfv,nstfv))

    dimtk=nomax*(nstfv-numin+1)
    allocate(mmatk(1:matsiz,1:dimtk))
    allocate(mmatkp(1:matsiz,1:dimtk))
    allocate(rmmat(1:matsiz,1:dimtk))
    
    allocate(tmat(1:matsiz,1:matsiz))
    allocate(m2tot(1:matsiz,1:matsiz))
    allocate(m2k(1:matsiz,1:matsiz))
    allocate(m2r(1:matsiz,1:matsiz))

    open(31,file='rotminm1.dat')
    write(31,*)'***** Mi_nm *****'
    open(32,file='rotminm2.dat')
    write(32,*)'***** Mi_nm^2 *****'

!-------------------------------------------------------------------
!   Direct BZ summation and check for M^i_{nm} rotation
!-------------------------------------------------------------------
    m2tot(:,:)=zzero
    do ik = 1, nkptnr
!
!     Read the matrix elements M^i_{nm}
!
      read(41,rec=ik) minmmat
      
      ijst=0
      do ist = 1, nomax
         do jst = numin, nstfv
            ijst=ijst+1
            mmatk(1:matsiz,ijst)=minmmat(1:matsiz,ist,jst)
          enddo ! jst
      enddo ! ist

      ikp=idikpq(indkpq(ik,iqp),iqp)
      read(41,rec=ikp) minmmat
      
      ijst=0
      do ist = 1, nomax
         do jst = numin, nstfv
            ijst=ijst+1
            mmatkp(1:matsiz,ijst)=minmmat(1:matsiz,ist,jst)
          enddo ! jst
      enddo ! ist
!
!     Get the matrix elements M^i_{nm} by symmetry
!
      isym=iksymq(ik,iqp)
!
!     Generate the MB rotation matrix according to the group symmetry operations
!
      call genmbrotmat(iq,isym)

!     Rotate M^i_{nm}
      call zgemm('c','n',matsiz,dimtk,matsiz, &
     &     zone,rotmat,matsiz,mmatkp,matsiz,zzero,rmmat,matsiz) 

      call boxmsg(31,'#','New K-point')
      write(31,*)
      write(31,*) 'ik = ', ik
      do jst = 1, ijst, 10
        do ist = 1, matsiz
          word="MT"; if(ist>locmatsiz)word="PW"
          write(31,'(a,2i5,6f12.6)') word, ist, jst, mmatk(ist,jst), rmmat(ist,jst), abs(mmatk(ist,jst)-rmmat(ist,jst))
        end do
      end do
      write(31,*)

!     calculate \sum_{nm} M^{i}_{nm} \times M^{j}_{nm}^{*} 
      call zgemm('n','c',matsiz,matsiz,dimtk, &
     &  zone,mmatk,matsiz,mmatk,matsiz,zzero,m2k,matsiz)
      
      m2tot(:,:)=m2tot(:,:)+m2k(:,:)

      call zgemm('n','c',matsiz,matsiz,dimtk, &
     &  zone,rmmat,matsiz,rmmat,matsiz,zzero,m2r,matsiz)

      call boxmsg(32,'#','New K-point')
      write(32,*)
      write(32,*) 'ik = ', ik
      do jst = 1, matsiz, 10
        do ist = 1, matsiz
          word="MT"; if(ist>locmatsiz)word="PW"
          write(32,'(a,2i5,6f12.6)') word, ist, jst, m2k(ist,jst), m2r(ist,jst), abs(m2k(ist,jst)-m2r(ist,jst))
        end do
      end do
      write(32,*)

    end do ! ik
    deallocate(mmatk)
    close(31)
    close(32)

!---------------------------------------------------------
!   \sum_{k}^{BZ} -> \sum_{k}^{BZ_q} \sum_{i}^{star(k)}
!---------------------------------------------------------

    open(36,file='summinm2.dat')
    write(36,*)'***** sum_k Mi_nm*Mj_nm *****'

    m2r(:,:)=zzero
    do ik = 1, nkptq(iqp)
      
      ikp = idikpq(ik,iqp)
      read(41,rec=ikp) minmmat
      
      ijst=0
      do ist = 1, nomax
         do jst = numin, nstfv
            ijst=ijst+1
            mmatkp(1:matsiz,ijst)=minmmat(1:matsiz,ist,jst)
          enddo ! jst
      enddo ! ist

      ! get M^i_{nm}*M^j_{nm}^{*}
      call zgemm('n','c',matsiz,matsiz,dimtk, &
     &     zone,mmatkp,matsiz,mmatkp,matsiz,zone,m2r,matsiz)
      
      ! sum over the symmetry operations which form the star of k
      do i = 2, nsymkstar(ik,iqp)
      
        isym = isymkstar(i,ik,iqp)

        ! Generate the MB rotation matrix according to the group symmetry operations
        call genmbrotmat(iq,isym)
        
        ! Rotate M^i_{nm}
        call zgemm('c','n',matsiz,dimtk,matsiz, &
       &     zone,rotmat,matsiz,mmatkp,matsiz,zzero,rmmat,matsiz) 

        ! M^i_{nm}*M^j_{nm}^{*}
        call zgemm('n','c',matsiz,matsiz,dimtk, &
       &     zone,rmmat,matsiz,rmmat,matsiz,zone,m2r,matsiz)
      
      end do ! i
    
    end do ! ik    

    call boxmsg(36,'#','Sum over k-points of Minm*Mjnm')
    write(36,*)
    do ist = 1, matsiz, 10
      do jst = 1, matsiz
        write(36,'(2i5,6f12.6)') ist, jst, m2tot(ist,jst), m2r(ist,jst), abs(m2tot(ist,jst)-m2r(ist,jst))
      end do
    end do
    close(36)

    deallocate(mmatkp)
    deallocate(rmmat) 
    deallocate(m2k,m2r,m2tot)
    
    deallocate(minmmat)
    close(41)

    if (input%gw%coreflag=='all') then
      
!--------------------------------------------
!   Check core-valence M^i_{cm} rotation
!--------------------------------------------

      inquire(IoLength=Recl) micmmat
      open(40,file='micmmat.io',action='READ',form='UNFORMATTED', &
   &       access='DIRECT',recl=Recl)

      if(allocated(micmmat))deallocate(micmmat)
      allocate(micmmat(locmatsiz,ncg,nstfv))
    
      dimtk=(nstfv-numin+1)*ncg
      allocate(mmatk(1:locmatsiz,1:dimtk))
      allocate(mmatkp(1:locmatsiz,1:dimtk))
      allocate(rmmat(1:locmatsiz,1:dimtk))
      allocate(m2k(1:locmatsiz,1:locmatsiz))
      allocate(m2r(1:locmatsiz,1:locmatsiz))
      mmatk=zzero
      mmatkp=zzero
      rmmat=zzero
      m2k=zzero
      m2r=zzero
      
      open(31,file='rotmicm1.dat')
      write(31,*)'***** Mi_cm *****'
      open(32,file='rotmicm2.dat')
      write(32,*)'***** Mi_cm^2 *****'
      
      do ik = 1, nkptnr
!
!       Read the matrix elements M^i_{cm}
!
        read(40,rec=ik) micmmat

        ijst=0
        do jst = numin, nstfv
          do icg = 1, ncg
            ijst=ijst+1
            mmatk(1:locmatsiz,ijst)=micmmat(1:locmatsiz,icg,jst)
          enddo ! icg
        enddo ! jst

        ikp=idikpq(indkpq(ik,iqp),iqp)
        read(40,rec=ikp) micmmat

        ijst=0
        do jst = numin, nstfv
          do icg = 1, ncg
            ijst=ijst+1
            mmatkp(1:locmatsiz,ijst)=micmmat(1:locmatsiz,icg,jst)
          enddo ! icg
        enddo ! jst
!
!       Get the matrix elements M^i_{cm} by symmetry
!
        isym=iksymq(ik,iqp)
!
!       Generate the MB rotation matrices according to the group symmetry operations
!
        call genmbrotmat(iq,isym)

!       Rotate M^i_{cm}
        call zgemm('n','n',locmatsiz,ijst,locmatsiz, &
       &     zone,rotmat(1:locmatsiz,1:locmatsiz),locmatsiz,mmatkp,locmatsiz,zzero,rmmat,locmatsiz)

        call boxmsg(31,'#','New K-point')
        write(31,*)
        write(31,*) 'ik = ', ik
        do jst = 1, ijst, 10
          do ist = 1, locmatsiz
            word="MT"
            write(31,'(a,2i5,6f12.6)') word, ist, jst, mmatk(ist,jst), rmmat(ist,jst), abs(mmatk(ist,jst)-rmmat(ist,jst))
          end do
        end do
        write(31,*)

!       calculate \sum_{nm} M^{i}_{nm} \times M^{j}_{nm}^{*} 
        call zgemm('n','c',locmatsiz,locmatsiz,ijst,zone,mmatk,locmatsiz,   &
       &            mmatk,locmatsiz,zzero,m2k,locmatsiz)
        
        call zgemm('n','c',locmatsiz,locmatsiz,ijst,zone,rmmat,locmatsiz,   &
       &            rmmat,locmatsiz,zzero,m2r,locmatsiz)

        call boxmsg(32,'#','New K-point')
        write(32,*)
        write(32,*) 'ik = ', ik
        do jst = 1, locmatsiz, 10
          do ist = 1, locmatsiz
            word="MT"
            write(32,'(a,2i5,6f12.6)') word, ist, jst, m2k(ist,jst), m2r(ist,jst), abs(m2k(ist,jst)-m2r(ist,jst))
          end do
        end do
        write(32,*)

      end do ! ik
      deallocate(mmatk)
      deallocate(mmatkp)
      deallocate(rmmat) 
      deallocate(m2k,m2r)
      
      close(31)
      close(32)
      
      deallocate(micmmat)
      close(40)
      
!--------------------------------------------
!   Check core-valence M^i_{nc} rotation
!--------------------------------------------

      inquire(IoLength=Recl) mincmat
      open(39,file='mincmat.io',action='READ',form='UNFORMATTED', &
   &       access='DIRECT',recl=Recl)

      if(allocated(mincmat))deallocate(mincmat)
      allocate(mincmat(locmatsiz,nstfv,ncg))
    
      dimtk=(nstfv-numin+1)*ncg
      allocate(mmatk(1:locmatsiz,1:dimtk))
      allocate(mmatkp(1:locmatsiz,1:dimtk))
      allocate(rmmat(1:locmatsiz,1:dimtk))
      allocate(m2k(1:locmatsiz,1:locmatsiz))
      allocate(m2r(1:locmatsiz,1:locmatsiz))
      mmatk=zzero
      mmatkp=zzero
      rmmat=zzero
      m2k=zzero
      m2r=zzero
      
      open(31,file='rotminc1.dat')
      write(31,*)'***** Mi_nc *****'
      open(32,file='rotminc2.dat')
      write(32,*)'***** Mi_nc^2 *****'
      
      do ik = 1, nkptnr
!
!       Read the matrix elements M^i_{nc}
!
        read(39,rec=ik) mincmat

        ijst=0
        do ist = numin, nstfv
          do icg = 1, ncg
            ijst=ijst+1
            mmatk(1:locmatsiz,ijst)=mincmat(1:locmatsiz,ist,icg)
          enddo ! icg
        enddo ! jst

        ikp=idikpq(indkpq(ik,iqp),iqp)
        read(39,rec=ikp) mincmat

        ijst=0
        do ist = numin, nstfv
          do icg = 1, ncg
            ijst=ijst+1
            mmatkp(1:locmatsiz,ijst)=mincmat(1:locmatsiz,ist,icg)
          enddo ! icg
        enddo ! jst
!
!       Get the matrix elements M^i_{cm} by symmetry
!
        isym=iksymq(ik,iqp)
!
!       Generate the MB rotation matrices according to the group symmetry operations
!
        call genmbrotmat(iq,isym)

!       Rotate M^i_{cm}
        call zgemm('n','n',locmatsiz,ijst,locmatsiz, &
       &     zone,rotmat(1:locmatsiz,1:locmatsiz),locmatsiz,mmatkp,locmatsiz,zzero,rmmat,locmatsiz)

        call boxmsg(31,'#','New K-point')
        write(31,*)
        write(31,*) 'ik = ', ik
        do jst = 1, ijst, 10
          do ist = 1, locmatsiz
            word="MT"
            write(31,'(a,2i5,6f12.6)') word, ist, jst, mmatk(ist,jst), rmmat(ist,jst), abs(mmatk(ist,jst)-rmmat(ist,jst))
          end do
        end do
        write(31,*)

!       calculate \sum_{nm} M^{i}_{nm} \times M^{j}_{nm}^{*} 
        call zgemm('n','c',locmatsiz,locmatsiz,ijst,zone,mmatk,locmatsiz,   &
       &            mmatk,locmatsiz,zzero,m2k,locmatsiz)
        
        call zgemm('n','c',locmatsiz,locmatsiz,ijst,zone,rmmat,locmatsiz,   &
       &            rmmat,locmatsiz,zzero,m2r,locmatsiz)

        call boxmsg(32,'#','New K-point')
        write(32,*)
        write(32,*) 'ik = ', ik
        do jst = 1, locmatsiz, 10
          do ist = 1, locmatsiz
            word="MT"
            write(32,'(a,2i5,6f12.6)') word, ist, jst, m2k(ist,jst), m2r(ist,jst), abs(m2k(ist,jst)-m2r(ist,jst))
          end do
        end do
        write(32,*)

      end do ! ik
      deallocate(mmatk)
      deallocate(mmatkp)
      deallocate(rmmat) 
      deallocate(m2k,m2r)
      
      close(31)
      close(32)

      deallocate(mincmat)
      close(39)
    
    end if ! core
    
    return
end subroutine
