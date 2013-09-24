!BOP
! !ROUTINE: init_relax
! !INTERFACE:
!
!
Subroutine init_relax
! !USES:
      Use modinput
      Use modmain
#ifdef TETRA
      Use modtetra
#endif
      Use modgw
      Use modxs
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, i1, i2, i3
      Integer :: ispn, nsym, isym, lspl
      Real (8) :: boxl (3, 4), vl (3), vc (3) 
      Real (8) :: blen(3), lambdab
      Integer(4), allocatable :: symmat(:,:,:)

!____________________________
! lattice and symmetry set up

      Call findsymcrys                      ! find the crystal symmetries and shift atomic positions if required

      Call findsymsite                      ! find the site symmetries

#ifdef XS
      Call findsymi &                       ! determine inverse symmery elements
     &  (input%structure%epslat, maxsymcrys, nsymcrys, symlat, lsplsymc, vtlsymc, isymlat, scimap)

      Call gensymt2 &                       ! generate symmetrization array for rank 2 tensors
     &  (maxsymcrys, nsymcrys, symlatc, lsplsymc, symt2)

      Call setupsym                         ! calculate advanced information on symmetry group
#endif

      Call checkmt                          ! check for overlapping muffin-tins

      If (allocated(sfacg)) &               ! allocate structure factor array for G-vectors
     &  deallocate (sfacg)
      Allocate (sfacg(ngvec, natmtot))

      Call gensfacgp &                      ! generate structure factors for G-vectors
     &  (ngvec, vgc, ngvec, sfacg)

      Call gencfun                          ! generate the characteristic function

      Call energynn                         ! determine the nuclear-nuclear energy

!_________________________________________________________________
! determine the k-point grid automatically from radkpt if required

      If (input%groundstate%autokpt) Then
          input%groundstate%ngridk(:) = Int &
         &  (input%groundstate%radkpt/&
         &  Sqrt(input%structure%crystal%basevect(1, :)**2+&
         &  input%structure%crystal%basevect(2, :)**2+&
         &  input%structure%crystal%basevect(3, :)**2)) + 1
      End If

!_____________________________________________________________________________
! if nktot is set (gt 0), determine the k-point grid automatically from nktot, 
! the total number of k-points

      If (input%groundstate%nktot.gt.0) Then
          blen(:)=sqrt(bvec(1,:)**2+bvec(2,:)**2+bvec(3,:)**2)           
          lambdab=Dble((input%groundstate%nktot/(blen(1)*blen(2)*blen(3)))**(1.d0/3.d0))
          input%groundstate%ngridk(:) = Max0(1, Int(lambdab*blen(:)+input%structure%epslat))
      End If

!______________________________
! setup the default k-point box

      boxl (:, 1) = input%groundstate%vkloff(:) / dble(input%groundstate%ngridk(:))
      boxl (:, 2) = boxl (:, 1)
      boxl (:, 3) = boxl (:, 1)
      boxl (:, 4) = boxl (:, 1)
      boxl (1, 2) = boxl (1, 2) + 1.d0
      boxl (2, 3) = boxl (2, 3) + 1.d0
      boxl (3, 4) = boxl (3, 4) + 1.d0

!________________________________________
! allocate the reduced k-point set arrays

      If (allocated(ivk)) deallocate (ivk)
      Allocate (ivk(3, input%groundstate%ngridk(1)*input%groundstate%ngridk(2)*input%groundstate%ngridk(3)))
      If (allocated(vkl)) deallocate (vkl)
      Allocate (vkl(3, input%groundstate%ngridk(1)*input%groundstate%ngridk(2)*input%groundstate%ngridk(3)))
      If (allocated(vkc)) deallocate (vkc)
      Allocate (vkc(3, input%groundstate%ngridk(1)*input%groundstate%ngridk(2)*input%groundstate%ngridk(3)))
      If (allocated(wkpt)) deallocate (wkpt)
      Allocate (wkpt(input%groundstate%ngridk(1)*input%groundstate%ngridk(2)*input%groundstate%ngridk(3)))
      If (allocated(ikmap)) deallocate (ikmap)
      Allocate (ikmap(0:input%groundstate%ngridk(1)-1, &
     &  0:input%groundstate%ngridk(2)-1, &
     &  0:input%groundstate%ngridk(3)-1))

!_________________________________
! generate the reduced k-point set

      if (input%groundstate%stypenumber < 0) then

!_____________________________________________________________
! suppress debug output in tetrahedron integration library (0)

          call tetrasetdbglv (0)

          nkpt = input%groundstate%ngridk(1)*input%groundstate%ngridk(2)*input%groundstate%ngridk(3)
          ntet = 6*nkpt
             
          if (allocated(indkp)) deallocate(indkp)
          allocate(indkp(nkpt))
          if (allocated(iwkp)) deallocate(iwkp)
          allocate(iwkp(nkpt))
          if (allocated(wtet)) deallocate(wtet)
          allocate(wtet(ntet))
          if (allocated(tnodes)) deallocate(tnodes)
          allocate(tnodes(4,ntet))

          nsym=1
          If (input%groundstate%reducek) nsym = nsymcrys

!__________________________________________             
! get rotational part of crystal symmetries

          if(allocated(symmat))deallocate(symmat)
          allocate(symmat(3,3,nsym))
          Do isym = 1, nsym
              lspl = lsplsymc(isym)

!_______________________________________________
! transpose of rotation for use with the library

              Do i1 = 1, 3
                  Do i2 = 1, 3
                      symmat(i1,i2,isym) = symlat(i2,i1,lspl)
                  End Do
              End Do
          End Do

          call factorize(3,input%groundstate%vkloff,ikloff,dkloff)

          call kgen(bvec,nsym,symmat,input%groundstate%ngridk,ikloff,dkloff,&
         &          nkpt,ivk,dvk,indkp,iwkp,ntet,tnodes,wtet,tvol,mnd)

!____________________
! getting ikmap array

          ik=0
          do i3=0,input%groundstate%ngridk(3)-1
          do i2=0,input%groundstate%ngridk(2)-1
          do i1=0,input%groundstate%ngridk(1)-1
              ik=ik+1
              ikmap(i1,i2,i3)=indkp(ik)
          end do
          end do
          end do

!_________________________________________________________              
! fractional and cartesian coordinates, and k-point weight

          do ik=1,nkpt
              vkl(:,ik)=dble(ivk(:,ik))/dble(dvk)
              call r3mv(bvec,vkl(:,ik),vkc(:,ik))
              wkpt(ik)=dble(iwkp(ik))/dble(input%groundstate%ngridk(1)* &
             &  input%groundstate%ngridk(2)*input%groundstate%ngridk(3))
          enddo ! ik
             
          deallocate(symmat)

      else

          Call genppts (input%groundstate%reducek, .False., &
         &              input%groundstate%ngridk, boxl, nkpt, ikmap, ivk, vkl, vkc, wkpt)

      end if ! tetra

!_____________________________
! G+k vectors: determine gkmax

      If ((input%groundstate%isgkmax .Ge. 1) .And. (input%groundstate%isgkmax .Le. nspecies)) Then
          gkmax = input%groundstate%rgkmax / rmt(input%groundstate%isgkmax)
      Else
          gkmax = input%groundstate%rgkmax / 2.d0
      End If
      If (2.d0*gkmax .Gt. input%groundstate%gmaxvr+input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(init1): 2*gkmax > gmaxvr  ", 2G18.10)') &
        & 2.d0 * gkmax, input%groundstate%gmaxvr
         Write (*,*)
         Stop
      End If

!_______________________________________
! find the maximum number of G+k-vectors

      Call getngkmax

!_______________________________
! allocate the G+k-vector arrays

      If (allocated(ngk)) deallocate (ngk)
      Allocate (ngk(nspnfv, nkpt))
      If (allocated(igkig)) deallocate (igkig)
      Allocate (igkig(ngkmax, nspnfv, nkpt))
      If (allocated(vgkl)) deallocate (vgkl)
      Allocate (vgkl(3, ngkmax, nspnfv, nkpt))
      If (allocated(vgkc)) deallocate (vgkc)
      Allocate (vgkc(3, ngkmax, nspnfv, nkpt))
      If (allocated(gkc)) deallocate (gkc)
      Allocate (gkc(ngkmax, nspnfv, nkpt))
      If (allocated(tpgkc)) deallocate (tpgkc)
      Allocate (tpgkc(2, ngkmax, nspnfv, nkpt))
      If (allocated(sfacgk)) deallocate (sfacgk)
      Allocate (sfacgk(ngkmax, natmtot, nspnfv, nkpt))

      Do ik = 1, nkpt
          Do ispn = 1, nspnfv
              vl (:) = vkl (:, ik)
              vc (:) = vkc (:, ik)

!_________________________
! generate the G+k-vectors

              Call gengpvec (vl, vc, ngk(ispn, ik), igkig(:, ispn, ik), vgkl(:, :, ispn, ik), &
             &               vgkc(:, :, ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, ik))

!___________________________________________
! generate structure factors for G+k-vectors

              Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), ngkmax, sfacgk(:, :, ispn, ik))

          End Do
      End Do

      Return
End Subroutine
!EOC
