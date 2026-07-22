module mod_get_eigenvectors_times_matchingcoefficients

contains

  !> Calculate the product of an eigenvector with the corresponding matching coefficients
  subroutine get_eigenvectors_times_matchingcoefficients(ik, trans, evec, evecalm)
    use modinput, only: input
    use modmain, only : ngkmax, apwordmax, lmmaxapw, natmtot, &
    &                   nspecies, natoms, idxas, idxlm, apword
    use modgw, only : Gkqset

    implicit none

    integer(4), intent(in) :: ik
    character(1), intent(in) :: trans
    complex(8), intent(in) :: evec(:,:)
    complex(8), intent(out) :: evecalm(:,:,:,:)

    integer(4) :: ia, is, ias
    integer(4) :: io
    integer(4) :: ist
    integer(4) :: l, m, lm
    integer(4) :: ngk
    complex(8), allocatable :: apwalm(:,:,:,:)

    complex(8), external :: zdotc
    complex(8), external :: zdotu

    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))

    ngk = Gkqset%ngk(1,ik)
    call match(ngk, &
    &          Gkqset%gkc(:,1,ik), &
    &          Gkqset%tpgkc(:,:,1,ik), &
    &          Gkqset%sfacgk(:,:,1,ik), &
    &          apwalm)

    select case (trans)
    case ('t', 'T')
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do l = 0, input%groundstate%lmaxapw
            do m = -l, l
              lm = idxlm(l,m)
              do io = 1, apword(l,is)
                do ist = 1, size(evec, 2)
                  evecalm(ist,io,lm,ias) = &
                  &  zdotu(ngk, &
                  &        evec(1:ngk,ist), 1, &
                  &        apwalm(1:ngk,io,lm,ias), 1)
                end do
              end do
            end do
          end do
        end do
      end do

    case ('c', 'C')
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do l = 0, input%groundstate%lmaxapw
            do m = -l, l
              lm = idxlm(l,m)
              do io = 1, apword(l,is)
                do ist = 1, size(evec, 2)
                  evecalm(ist,io,lm,ias) = &
                  &  zdotc(ngk, &
                  &        apwalm(1:ngk,io,lm,ias), 1, &
                  &        evec(1:ngk,ist), 1)
                end do
              end do
            end do
          end do
        end do
      end do

    case default
      write(*,*) 'ERROR in get_eigenvectors_times_matchingcoefficients'
      write(*,*) 'Allowed values of trans are t, T, c or C'
      write(*,'(a19,a1)') 'Received trans = ', trans
      stop
    end select

    deallocate(apwalm)
  end subroutine get_eigenvectors_times_matchingcoefficients

end module mod_get_eigenvectors_times_matchingcoefficients
