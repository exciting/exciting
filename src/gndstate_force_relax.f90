subroutine gndstate_force_relax(redoscl)
  use modmain
  use modmpi
  logical ,intent(out)::redoscl
  if (rank.eq.0)then
     redoscl=.false.
     write(60,*)
     write(60,'("+------------------------------+")')
     write(60,'("| Self-consistent loop stopped |")')
     write(60,'("+------------------------------+")')
     ! write density and potentials to file only if maxscl > 1
     if (maxscl.gt.1) call writestate
  endif
  !-----------------------!
  !     compute forces    !
  !-----------------------!
  if ((.not.tstop).and.(tforce)) then
     call force
     if (rank.eq.0) then
        ! output forces to INFO.OUT
        call writeforce(60)
        ! write maximum force magnitude to FORCEMAX.OUT
        write(64,'(G18.10)') forcemax
        call flushifc(64)
     endif
  end if
  !---------------------------------------!
  !     perform structural relaxation     !
  !---------------------------------------!
  if ((.not.tstop).and.((task.eq.2).or.(task.eq.3))) then
     if (rank.eq.0) then
        write(60,*)
        write(60,'("Maximum force magnitude (target) : ",G18.10," (",G18.10,")")') &
             forcemax,epsforce
        call flushifc(60)
     endif
     ! check force convergence
     if (forcemax.le.epsforce) then
        if (rank.eq.0) then
           write(60,*)
           write(60,'("Force convergence target achieved")')
        endif
        goto 30
     end if
     ! update the atomic positions if forces are not converged
     call updatpos
     if (rank.eq.0) then
        write(60,*)
        write(60,'("+--------------------------+")')
        write(60,'("| Updated atomic positions |")')
        write(60,'("+--------------------------+")')
        do is=1,nspecies
           write(60,*)
           write(60,'("Species : ",I4,", ",A)') is,trim(spsymb(is))
           write(60,'(" atomic positions (lattice) :")')
           do ia=1,natoms(is)
              write(60,'(I4,3F14.8)') ia,atposl(:,ia,is)
           end do
        end do
        ! write lattice vectors and optimised atomic positions to file
        call writegeom(.true.)
     endif
     ! check for overlapping muffin-tins
     call checkmt
     ! generate structure factors for G-vectors
     call gensfacgp(ngvec,vgc,ngvec,sfacg)
     ! generate the characteristic function
     call gencfun
     ! generate structure factors for G+k-vectors
     do ispn=1,nspnfv
        do ik=1,nkpt
           call gensfacgp(ngk(ik,ispn),vgkc(1,1,ik,ispn),ngkmax,sfacgk(1,1,ik,ispn))
        end do
     end do
     ! determine the new nuclear-nuclear energy
     call energynn
     ! add blank line to TOTENERGY.OUT, FERMIDOS.OUT and MOMENT.OUT
     if (rank.eq.0) then
        write(61,*)
        write(62,*)
        if (spinpol) write (63,*)
     endif
     ! begin new self-consistent loop with updated positions
     redoscl=.true.

  end if

30 continue
  !bcasts!!

#ifdef MPI
  call MPI_bcast(redoscl,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(fsum,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(atposc,size(atposc),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(atposl,size(atposl),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vkl,size(vkl),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vkc,size(vkc),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(sfacg,size(sfacg),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngk,size(ngk),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vgkc,size(vgkc),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(sfacgk,size(sfacgk),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngkmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(engynn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

#endif

end subroutine gndstate_force_relax
