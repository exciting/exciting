subroutine borncharges( pol, disp, zstar, sumrule, mode)
  use modmpi
  use modinput
  use mod_lattice
  use mod_atoms

  real(8), intent( in)      :: pol(3,2,3,natmtot), disp
  logical, intent( in)      :: sumrule
  character(*), intent( in) :: mode
  real(8), intent( inout)   :: zstar(3,3,0:natmtot)

  integer :: ia, is, ias, ip, vi(3), un
  real(8) :: vr(3)

  if( trim( mode) == 'write') then
    call put( zstar)
    return
  elseif( trim( mode) == 'read') then
    call get( zstar)
    return
  end if

  zstar = 0.d0
  do is = 1, nspecies
    do ia = 1, natoms( is)
      ias = idxas( ia, is)
      do ip = 1, 3
        call r3mv( ainv, pol(:,2,ip, ias)-pol(:,1,ip,ias), vr)
        vr = vr*omega
        call r3ws( 1.d-16, input%structure%crystal%basevect, vr, vi)
        call r3mv( input%structure%crystal%basevect, vr, zstar(ip,:,ias))
      end do
      zstar(:,:,ias) = zstar(:,:,ias)/disp
      zstar(:,:,0) = zstar(:,:,0) + zstar(:,:,ias)
    end do
  end do
  zstar(:,:,0) = zstar(:,:,0)/natmtot

  if( sumrule) then
    do ias = 1, natmtot
      zstar(:,:,ias) = zstar(:,:,ias) - zstar(:,:,0)
    end do
  end if
  
  return

  contains
    subroutine put( z)
      use mod_misc, only: filext
      use m_getunit
      real(8), intent( in) :: z(3,3,0:natmtot)

      integer :: un

      if( mpiglobal%rank == 0) then
        call getunit( un)
        open( un, file='ZSTAR'//trim( filext), action='write', form='formatted')
        write( un, '("# Born effective charges")')
        if( sumrule) write( un, '("# Note: Acoustic sum rule was automatically imposed.")')
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            call r3mv( ainv, atposc(:,ia,is), vr)
            write( un,'("# Species ",i2,"  Atom ",i2," (",a2,i2") : ",3f13.6)') is, ia, spsymb( is), ia, vr
            do ip = 1, 3
              write( un, '(3f16.6)') z(ip,:,ias)
            end do
          end do
        end do
        if( sumrule) then
          write( un, '("# acoustic sum rule correction")')
          do ip = 1, 3
            write( un, '(3f16.6)') z(ip,:,0)
          end do
        end if
        close( un)
      end if
      call barrier

      return
    end subroutine

    subroutine get( z)
      use m_getunit
      use mod_misc, only: filext
      real(8), intent( out) :: z(3,3,0:natmtot)
      
      integer :: un, i, l, ias
      character(256) :: fname, buf
      logical :: exist

      z = 0.d0
      write( fname, '("ZSTAR")')
      fname = trim( fname)//trim( filext)
      inquire( file=trim( fname), exist=exist)
      if( .not. exist) then
        if( mpiglobal%rank == 0) then
          write(*,*)
          write(*,'("Error (borncharges): File ",a," does not exist.")') trim( fname)
        end if
        call terminate
      end if
      call getunit( un)
      open( un, file=trim( fname), status='old', form='formatted')
      l = 0; ias = 0
      do
        read( un, '(a)', iostat=i) buf
        if( i /= 0) exit
        if( (buf( 1:1) == '#') .or. (len( trim( buf)) == 0)) cycle
        if( modulo( l, 3) == 0) ias = ias + 1
        if( ias <= natmtot) then
          i = modulo( l, 3) + 1
          read( buf, *) z(i,1,ias), z(i,2,ias), z(i,3,ias)
        end if
        if( ias == natmtot+1) then
          i = modulo( l, 3) + 1
          read( buf, *) z(i,1,0), z(i,2,0), z(i,3,0)
        end if
        l = l + 1
      end do
      close( un)

      if( (ias > natmtot) .and. .not. sumrule) then
        do ias = 1, natmtot
          z(:,:,ias) = z(:,:,ias) + z(:,:,0)
        end do
      end if
      if( sumrule) ias = ias - 1
      if( ias < natmtot) then
        if( mpiglobal%rank == 0) then
          write(*,*)
          write(*,'("Error (borncharges): File ",a," contains not enough data.")') trim( fname)
        end if
        call terminate
      end if

      return
    end subroutine
end subroutine borncharges
