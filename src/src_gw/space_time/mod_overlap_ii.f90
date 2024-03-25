module mod_overlap_ii
    use precision, only: wp
    use modinput, only: input
    use mod_lattice, only: omega
    use modgw, only: Gqset, Gset
    use mod_gvector, only: cfunig 
    use mod_product_basis, only: sgi
    use constants, only: zzero, zone
    ! !! delete if not needed
    ! use modgw, only: Gkqset
    use modgw, only: mpwipw

    implicit none
    private
    public :: overlap_ii

    contains

    !!--------------------------------------------------------------------------------------------------
    !! Overlap matrix in MT-I region 
    !> \begin{equation}
    !>   \left\langle  e^{i(\mathbf{q}+ \mathbf{G}^{\prime})\mathbf{r}^{\prime} } 
    !>   \mid \chi^{\mathbf{q}}_{I^{\prime}}\right\rangle_{IR} = \left\langle  
    !>   e^{i(\mathbf{q}+ \mathbf{G}^{\prime})\mathbf{r}^{\prime} } \mid 
    !>   e^{i(\mathbf{q}+ \mathbf{G}_{I^{\prime}})\mathbf{r}^{\prime} } \right\rangle_{IR} = 
    !>   \int d\mathbf{r}^{\prime} e^{-i(\mathbf{G}' - \mathbf{G}_{I^{\prime}})\mathbf{r}^{\prime} }  
    !>   = \Theta(\mathbf{G}' - \mathbf{G}_{I^{\prime}}) = S_{\mathbf{G}',I'}^{\mathbf{q}}
    !> \end{equation}
    !! TODO: Is above equation is correct or do I need eqn(40 of FHI-gap paper in place of $chi ???
    !!--------------------------------------------------------------------------------------------------

    subroutine overlap_ii(iq, s_ii)
        !> q-point index
        integer, intent(in) :: iq
        !> Overlap matrix: I-I case
        complex(wp), allocatable, intent(out) :: s_ii(:,:)

        !> Local variables 
        integer :: ig, ig1, ig2, ngq
        integer :: igv(3) !! integer coodintates of G-G'
        !> Step function matrix in interstitial region \[ \Theta(G-G') \]
        complex(wp), allocatable :: step_fn_mat(:,:)


        ngq = Gqset%ngk(1, iq)   !!! TODO: do I need Gkqset or Gqset ???
        
        if(allocated(s_ii)) deallocate(s_ii)
        allocate(s_ii(ngq, ngq))
        allocate(step_fn_mat(ngq, ngq))
        s_ii(:,:) = zzero
        step_fn_mat(:,:) = zzero

        call gencfun   !! This gives heaviside step function in "cfunig(:)" array
    
        !! TODO: Do I need omega(volume) or it's taken care of in "call gencfun" ?????????
        !! TODO: Do I need "conjg(cfunig(:))" or only "cfunig(:)"? (sgi* or sgi of FHI-gap)

        ! Old
        ! do ig1 = 1, ngq
        !   ig = Gset%ivgig(0,0,0)
        !   s_ii(ig1, ig1) = conjg(cfunig(ig))   !! diagonal terms
        !   do ig2 = ig1+1, ngq
        !     igv(:) = Gset%ivg(:, Gqset%igkig(ig1,1,iq)) - Gset%ivg(:, Gqset%igkig(ig2,1,iq))
        !     ig = Gset%ivgig(igv(1),igv(2),igv(3))
        !     s_ii(ig2, ig1) = cfunig(ig)           !! off-diagonal terms (lower triangle)
        !     s_ii(ig1, ig2) = conjg(cfunig(ig))    !! off-diagonal terms (upper triangle)
        !   enddo
        ! enddo

        do ig1 = 1, ngq
          ig = Gset%ivgig(0,0,0)
          step_fn_mat(ig1, ig1) = conjg(cfunig(ig))   !! diagonal terms
          do ig2 = ig1+1, ngq
            igv(:) = Gset%ivg(:, Gqset%igkig(ig1,1,iq)) - Gset%ivg(:, Gqset%igkig(ig2,1,iq))
            ig = Gset%ivgig(igv(1),igv(2),igv(3))
            step_fn_mat(ig2, ig1) = cfunig(ig)           !! off-diagonal terms (lower triangle)
            step_fn_mat(ig1, ig2) = conjg(cfunig(ig))    !! off-diagonal terms (upper triangle)
          enddo
        enddo

      call diagsgi(iq)

      ! call calcmpwipw(iq)   !! TEST and DELETE
      ! do ig1 = 1, size(mpwipw, 1)
      !   write(16,*) ig1, real(mpwipw(ig1, 1)), imag(mpwipw(ig1, 1))
      ! enddo


      ! call zgemm('n', 'n', ngq, ngq, ngq,  zone, step_fn_mat, ngq,  sgi, ngq,  zzero, s_ii, ngq)
      call zgemm('n', 'n', ngq, ngq, ngq,  dsqrt(omega), step_fn_mat, ngq,  sgi, ngq,  zzero, s_ii, ngq)

       do ig1 = 1, ngq
         write(16,*) ig1, real(s_ii(ig1, 1)), real(sgi(ig1, 1)), real(step_fn_mat(ig1, 1))!, imag(s_ii(ig1, 1))
         write(17,*) ig1, real(s_ii(ig1, 3)), real(sgi(ig1, 3)), real(step_fn_mat(ig1, 3))!, imag(s_ii(ig1, 1))
         write(18,*) ig1, real(s_ii(1, ig1)), real(sgi(1, ig1)), real(step_fn_mat(1, ig1))!, imag(s_ii(ig1, 1))
         write(19,*) ig1, real(cfunig(ig1)), real(sgi(ig1, 1)), real(step_fn_mat(ig1, 1))!, imag(cfunig(ig1))
       enddo

    end subroutine overlap_ii

end module mod_overlap_ii 
