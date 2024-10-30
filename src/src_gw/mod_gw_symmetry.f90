!> This module contains the subroutines and procedures related to symmetry operations in the mixed basis and the Coulomb eigenbasis.
!> The overall idea is that of Eq. 51 of https://doi.org/10.1016/j.cpc.2011.12.006. In a nutshell, 
!> the inverse of the dielectric matrix (or any other quantity in the mixed basis or the Coulomb 
!> eigenbasis) at a general q-point is recovered from its representative (i.e., irreducible q) by applying a matrix multiplication representing the symmetry operation
!> that relates the representative to the desired q-point.
module mod_gw_symmetry

    use precision, only: dp, i32
    use mod_kpointset, only: k_set, kq_set, G_set, Gk_set

    implicit none

    private
    public :: rotate_from_qa_to_qb_matrices

contains

    !> This subroutine computes the (adjoint) rotation matrix in the mixed basis, $\mathbf{R}$, mapping the point
    !> iqa (\[\mathbf{q}_{a}\]) to iqb (\[\mathbf{q}_{b}\]). Therefore, given a vector at iqa expressed in terms of the mixed basis,
    !> i.e. \[\mathbf{v}(\mathbf{q}_{a})\], we can obtain the symmetry related vector for iqb by 
    !> \[\mathbf{v}(\mathbf{q}_{b}) = \mathbf{R}^\dagger \mathbf{v}(\mathbf{q}_{a})\]\]. 
    subroutine compute_rotation_matrix_mb(iqa, iqb, kset, kqset, Gqset, Gset, gw_io_format, rotmat)

        use constants,         only: zzero
        use mod_lattice,       only: bvec, binv
        use mod_symmetry,      only: nsymcrys, symlat, lsplsymc, vtlsymc, symlatc, what_maps_q_2_qp
        use mod_product_basis, only: matsiz, locmatsiz

        implicit none

        !> The index to qa
        integer(i32), intent(in)                :: iqa
        !> The index to qb
        integer(i32), intent(in)                :: iqb
        !> Reciprocal mesh
        type(k_set), intent(in)                 :: kset
        !> Reciprocal double mesh
        type(kq_set), intent(in)                :: kqset
        !> G + k but with Gmax from GW not GS
        type(Gk_set), intent(in)                :: Gqset
        !> G-mesh set
        type(G_set),  intent(in)                :: Gset
        !> The format in which the files format has been saved
        character(len=*), intent(in)            :: gw_io_format
        !> Rotation matrix in MB
        complex(dp), allocatable, intent(out) :: rotmat(:,:)

        ! Symmetry operation mapping qa into qb
        integer(i32) :: isym

        ! Symmetry operation in reduced and cartesian coordinates (rotation)
        integer(i32) :: rot_reduced(3,3) 
        real(dp)     :: rot_cartesian(3,3)
        ! The translation in reduced coordinates
        real(dp) :: tau_reduced(3)

        ! integer to contain the symmetry operation id in full list
        integer(i32) :: lspl

        ! G-vector needed to map the q into the first BZ
        real(dp) :: G_red(3)

        ! Get the proper size for the rotmat
        matsiz = max(Gqset%ngk(1,iqa),Gqset%ngk(1,iqb)) + locmatsiz

        ! Allocate rotmat (intent(out) always deallocates an array if it is allocatable argument)
        allocate(rotmat(matsiz, matsiz), source=zzero)

        ! Get the symmetry operation that maps qb into qa
        isym = what_maps_q_2_qp(iqb, iqa, kqset%vql, nsymcrys, lsplsymc, symlatc, bvec, binv, &
                                kset%ngridk, kset%vkloff)

        ! Notice we are looking for the inverse of the symmetry operation
        lspl               = lsplsymc(isym)
        rot_reduced(:,:)   = symlat(:,:,lspl)
        rot_cartesian(:,:) = symlatc(:,:,lspl)
        tau_reduced(:)     = vtlsymc(:,isym)

        ! Get the G vector such that after the rotation is applied to iqb it is mapped into 
        ! the first BZ. Notice, that this definition differs from https://doi.org/10.1016/j.cpc.2011.12.006 in the fact
        ! that we are using the inverse operation and a minus. 
        ! \mathbf{q}_a = \mathbf{R}^{a\rightarrow b} \mathbf{q}_b - \mathbf{G}
        ! Moreover, the expression differs with the top one as exciting rotations are different
        ! from the spglib description (i.e. translation is applied before rotation and by the opposite side).  
        G_red(:)  = matmul(kqset%vql(:,iqb), real(rot_reduced, kind=dp)) - kqset%vql(:,iqa)

        ! Computes the MT contribution to the MB representation of the symmetry operation
        call compute_rotation_matrix_mb_mt(iqa, iqb, kqset, isym, rot_cartesian, G_red, tau_reduced, rotmat)

        ! Computes the IPW contribution to the MB representation of the symmetry operation
        call compute_rotation_matrix_mb_ipw(iqa, iqb, Gqset, Gset, gw_io_format, rot_reduced, G_red, tau_reduced, rotmat)
        
    end subroutine compute_rotation_matrix_mb

    !> Computes the MT contribution to the MB representation of the symmetry operation
    subroutine compute_rotation_matrix_mb_mt(iqa, iqb, kqset, isym, rot_cartesian, G_red, tau_reduced, rotmat)

        use constants,          only: twopi
        use mod_symmetry,       only: ieqatom
        use mod_misc_gw,        only: atposl
        use mod_atoms,          only: idxas
        use mod_product_basis,  only: locmatsiz, mbindex

        implicit none

        !> The index to qa
        integer(i32), intent(in)                :: iqa
        !> The index to qb
        integer(i32), intent(in)                :: iqb
        !> Reciprocal double mesh
        type(kq_set), intent(in)                :: kqset
        !> Index of the symmetry operation that maps qb into qa
        integer(i32), intent(in)                :: isym
        !> Symmetry operation in reduced and cartesian coordinates 
        !> Symmetry operation (rotation) in reduced coordinates
        real(dp), intent(in)                    :: rot_cartesian(3,3) 
        !> Get the G vector such that after the rotation is applied to iqb from iqa, the former it is mapped into 
        !> the first BZ.
        real(dp), intent(in)                    :: G_red(3)
        !> Symmetry operation (translation) in reduced coordinates
        real(dp), intent(in)                    :: tau_reduced(3)
        !> The rotation matrix to fill
        complex(dp), intent(inout)              :: rotmat(:,:)
        

        ! Integer indexes 
        integer(i32) :: imix, jmix
        ! Indexes for MT part
        integer(i32) :: imix_is, imix_ia, imix_irm, imix_bl, imix_bm
        integer(i32) :: ja, jas
        integer(i32) :: jmix_is, jmix_ia, jmix_irm, jmix_bl, jmix_bm
        ! Product of G vector with an atomic position
        real(dp) :: Gatpos
        ! Product of a G vector and a translation
        real(dp) :: Gtau
        
        ! Rotation to the spherical harmonic
        complex(dp), external :: getdlmm

        ! \[ R_{ji} = e^{2\pi\mathbf{G^{q}_{R}}\cdot\mathbf{r}_{\beta}} e^{-2\pi\mathbf{{q}_{b}}\cdot\mathbf{\tau}} D^L_{m_i m_j} \]
        ! The first term comes from the rotation of \[e^{i\mathbf{q}_b \cdot r_\alpha}\], the second comes from a 
        ! general phase arising in the rotation of the Bloch factor (see https://docs.abinit.org/theory/wavefunctions/), and the last term is the Wigner D-matrix coming from the rotation
        ! of the spherical harmonics.

        ! We compute it using threads, notice that due to branching the computations
        ! are unbalanced, thus the dynamic execution. 

        !$omp parallel do default(none) private(imix, imix_is, imix_ia, imix_irm, imix_bl, imix_bm) &
        !$omp private(jmix, jmix_is, jmix_ia, jmix_irm, jmix_bl, jmix_bm, ja, jas) &
        !$omp private(Gatpos) shared(locmatsiz, mbindex, ieqatom, idxas, G_red, atposl, Gtau) &
        !$omp shared(rot_cartesian, tau_reduced, rotmat, kqset, isym, iqb) schedule(dynamic)
        do imix = 1, locmatsiz
            ! Get information about the mixed basis
            imix_is  = mbindex(imix,1)
            imix_ia  = mbindex(imix,2)
            imix_irm = mbindex(imix,3)
            imix_bl  = mbindex(imix,4)
            imix_bm  = mbindex(imix,5)

            ! Get the equivalent atom due to the symmetry operation
            ja  = ieqatom(imix_ia,imix_is,isym)
            jas = idxas(ja,imix_is)

            ! Iterate over the other MT mixed basis componenents
            do jmix = 1, locmatsiz  
                jmix_is  = mbindex(jmix,1)
                jmix_ia  = mbindex(jmix,2)
                jmix_irm = mbindex(jmix,3)
                jmix_bl  = mbindex(jmix,4)
                jmix_bm  = mbindex(jmix,5)

                ! If not the atom we are interested in skip 
                if (jmix_is /= imix_is .or. ja /= jmix_ia) cycle
                
                ! If not the same radial function and bigl skip
                if (jmix_irm /= imix_irm .or. jmix_bl /= imix_bl) cycle

                ! Compute the term coming from the rotation of the exponential 
                ! factor with the atomic position (i.e. \[e^{i\mathbf{q}_a \cdot r_\beta}\]). Notice 
                ! the the possible lattice vector mapping to the equivalent \[ r_\alpha\], would equal to 
                ! 1 as G is a reciprocal lattice vector.
                Gatpos =  twopi * dot_product(G_red(:), atposl(:,ja,imix_is))

                ! Compute the rotation term coming from rotation of Bloch like terms
                Gtau   = -twopi * dot_product(kqset%vql(:,iqb), tau_reduced)

                ! Compute the contribution coming from the rotation of the spherical harmonics 
                ! getdlmm computes the Wigner-D matrix element
                ! See getdlmm.f90 for further information regarding that term.  
                rotmat(jmix,imix) = exp(cmplx(0.0_dp, Gatpos, kind=dp)) * exp(cmplx(0.0_dp, Gtau, kind=dp)) * &
                                    getdlmm(rot_cartesian, imix_bl, jmix_bm, imix_bm)

            end do

        end do
        !$omp end parallel do

    end subroutine compute_rotation_matrix_mb_mt

    !> Computes the IPW contribution to the MB representation of the symmetry operation
    subroutine compute_rotation_matrix_mb_ipw(iqa, iqb, Gqset, Gset, gw_io_format, rot_reduced, G_red, tau_reduced, rotmat)

        use constants,                     only: twopi
        use mod_Gvector,                   only: cfunig
        use mod_product_basis,             only: locmatsiz, basename_sgi
        use gw_io,                         only: read_from_file, build_file_name
        use general_matrix_multiplication, only: matrix_multiply

        implicit none

        !> The index to qa
        integer(i32), intent(in)                :: iqa
        !> The index to qb
        integer(i32), intent(in)                :: iqb
        !> G + k but with Gmax from GW not GS
        type(Gk_set), intent(in)                :: Gqset
        !> G-mesh set
        type(G_set),  intent(in)                :: Gset
        !> The format in which the files format has been saved
        character(len=*), intent(in)            :: gw_io_format
        !> Symmetry operation (rotation) in reduced coordinates
        integer(i32), intent(in) :: rot_reduced(3,3) 
        !> Get the G vector such that after the rotation is applied to iqb from iqa, the former it is mapped into 
        !> the first BZ.
        real(dp), intent(in)     :: G_red(3)
        !> Symmetry operation (translation) in reduced coordinates
        real(dp), intent(in)     :: tau_reduced(3)
        !> The rotation matrix to fill
        complex(dp), intent(inout) :: rotmat(:,:)

        ! Integer indexes 
        integer(i32) :: imix, jmix
        ! Indexes for IPW part
        integer(i32) :: ngq_a, ngq_b, igqb, igqa, ig
        !  G-vectors
        integer(i32) :: g1(3), g2(3), g3(3)
        ! Filename for the sgi
        character(len=30)   :: sgi_file_name_a, sgi_file_name_b
        ! Sgi components for each q-point
        complex(dp), allocatable :: sgi_a(:,:), sgi_b(:,:)
        ! Temporaries for the IPW part
        complex(dp), allocatable :: tmat1(:,:), tmat2(:,:)
        ! Product of a G vector and a translation
        real(dp) :: Gtau

        ! \[ R_{ji} = \sum_{G'} \mathbf{S}^{*}_{\mathf{G},j}(\mathbf{q}_b)  \sum_{G} e^{-i(\mathbf{G}+\mathbf{q}_b)\tau} 
        ! \mathcl{l}_{\mathbf{RG-G'+G_{R}} \mathbf{S}_{\mathf{G},i}(\mathbf{q}_a) \]
        !
        ! Notice this equivalent to the usual PW terms (see Eq 51 in https://doi.org/10.1016/j.cpc.2011.12.006) 
        ! but including the overlap between PW comming from the specifics of our basis set.
        ngq_a = Gqset%ngk(1,iqa)
        ngq_b = Gqset%ngk(1,iqb)
        
        allocate(tmat1(ngq_b, ngq_a))
        allocate(tmat2(ngq_b, ngq_b))

        !$omp parallel do collapse(2) default(none) private(igqb, igqa, g1, g2, g3, ig, Gtau) &
        !$omp shared(ngq_a, ngq_b, Gset, Gqset, iqa, iqb, rot_reduced, G_red, tau_reduced, cfunig, tmat1)
        do igqa = 1, ngq_a  ! loop over qa+G
            do igqb = 1, ngq_b   ! loop over qb+G'

                g1(:) = Gset%ivg(:,Gqset%igkig(igqb,1,iqb))
                g2(:) = Gset%ivg(:,Gqset%igkig(igqa,1,iqa))
                g3(:) = matmul(g1(:),rot_reduced) - g2(:) + nint(G_red(:))
                ! Index of (RG-G'+G_{R}) vector
                ig = Gset%ivgig(g3(1),g3(2),g3(3)) 

                Gtau = -twopi * dot_product(Gqset%vgkl(:,igqb,1,iqb),tau_reduced)

                ! The first is the usual expression for PW, the second comes from the fact 
                ! that we are not working with PW but IPW, so we need to add the Heaviside function
                ! which zeroes inside the MT.
                tmat1(igqb,igqa) = exp(cmplx(0.0_dp, Gtau, kind=dp)) * conjg(cfunig(ig))

            end do 
        end do
        !$omp end parallel do
        
        ! Here load the Sgi for each q point    
        call build_file_name( basename_sgi, iqa, sgi_file_name_a )
        call read_from_file( sgi_file_name_a, sgi_a, gw_io_format )
        call build_file_name( basename_sgi, iqb, sgi_file_name_b )
        call read_from_file( sgi_file_name_b, sgi_b, gw_io_format )

        ! Including the Sgi from iqa and iqb into the rotation terms accounting
        ! for the IPW terms of the mixed basis
        call matrix_multiply(tmat1, sgi_a, tmat2,'n','n')
        
        tmat1 = reshape(tmat1, [ngq_a, ngq_b])
        ! This constructs tells the computer to either use SIMD and/or threads to 
        ! perform this unstructured operation
        !$omp parallel workshare 
        sgi_b(:,:) = conjg(sgi_b(:,:))
        !$omp end parallel workshare
        call matrix_multiply(tmat2, sgi_b, tmat1,'t','n')

        ! Force the use of SIMD instructions for copying
        ! if unsupported it falls back to normal threads and simd is ignored
        !$omp parallel do simd collapse(2) 
        do imix = 1, ngq_a
            do jmix = 1, ngq_b
              rotmat(locmatsiz+jmix,locmatsiz+imix) = tmat1(jmix,imix)
            end do
        end do
        !$omp end parallel do simd 

    end subroutine compute_rotation_matrix_mb_ipw

    !> This subroutine changes the basis for a rotation matrix from the mixed basis
    !> to the appropriate Coulomb eigenbasis. Notice this will fail if iqa .and. iqb are Gamma, as we are not 
    !> removing the eigenvector corresponding to a plane-wave with \(G=0\) in the case of \(q=0\).
    subroutine rotmatmb_2_bare_coulomb_eigenbasis(iqa, iqb, kqset, barcevtol, gw_io_format, rotmat_mb, rotmat_cb)

        use mod_coulomb_potential, only: vmat, barcev, read_barcev_vmat_from_file, &
                                         index_of_first_element_above_threshold
        use mod_product_basis,     only: mbsiz, matsiz
        use general_matrix_multiplication, only: matrix_multiply
        use modmpi,                only: terminate
        use mod_misc_gw,           only: gammapoint

        implicit none 

        !> The index to qa
        integer(i32), intent(in)                :: iqa
        !> The index to qb
        integer(i32), intent(in)                :: iqb
        !> Reciprocal double mesh
        type(kq_set), intent(in)                :: kqset
        !> The tolerance factor is used to reduce the size of the Vc-diagonal product basis
        real(dp),     intent(in)                :: barcevtol
        !> The format in which the files format has been saved
        character(len=*), intent(in)            :: gw_io_format  
        !> Rotation matrix in mixed basis
        complex(dp), intent(in)                 :: rotmat_mb(:,:) 
        !> Rotation matrix in Coulomb eigenbasis
        complex(dp), allocatable, intent(out)   :: rotmat_cb(:,:)

        ! Coulomb eigenvectors
        complex(dp), allocatable :: barc_eigenvectors_qa(:,:), barc_eigenvectors_qb(:,:)

        ! index
        integer(i32) :: i, idx

        complex(dp), allocatable :: temporary(:,:)

        ! As index_of_first_element_above_threshold does not remove eigenvector corresponding to a plane-wave with \(G=0\) 
        ! in the case of \(q=0\), we call terminate if working with the Gamma point
        if (gammapoint(kqset%vqc(:,iqa),tol=1.0e-6_dp) .or. gammapoint(kqset%vqc(:,iqb),tol=1.0e-6_dp)) then
            call terminate('Error(rotmatmb_2_bare_coulomb_eigenbasis): This function is not suitable for the Gamma point')
        end if

        call read_barcev_vmat_from_file(iqa, gw_io_format)
        idx = index_of_first_element_above_threshold(barcev, barcevtol)
        allocate(barc_eigenvectors_qa, source=vmat(:,idx:))
        
        call read_barcev_vmat_from_file(iqb, gw_io_format)
        idx = index_of_first_element_above_threshold(barcev, barcevtol)
        allocate(barc_eigenvectors_qb, source=vmat(:,idx:))

        ! The read overwrites the correct values
        matsiz = size(barc_eigenvectors_qb, 1)
        mbsiz  = size(barc_eigenvectors_qb, 2)

        ! Transforming the rotation matrix to the bare Coulomb eigenbasis
        allocate(temporary(matsiz,mbsiz))
        allocate(rotmat_cb(mbsiz,mbsiz))
        
        call matrix_multiply(rotmat_mb, barc_eigenvectors_qb, temporary)
        call matrix_multiply(barc_eigenvectors_qa, temporary, rotmat_cb, 'c', 'n')

    end subroutine rotmatmb_2_bare_coulomb_eigenbasis

    !> This subroutine maps the matrix in Coulomb eigenbasis from qa to qb
    !> using crystal symmetry. It raises an error in case no symmetry operation
    !> maps qa to qb.
    subroutine rotate_from_qa_to_qb_matrices(iqa, iqb, barcevtol, gw_io_format, matrix_qa, matrix_qb)

        use modgw, only: kset, kqset, Gqset, Gset, time_rotmb
        use general_matrix_multiplication, only: matrix_multiply

        implicit none

        !> The index to qa
        integer(i32), intent(in)                :: iqa
        !> The index to qb
        integer(i32), intent(in)                :: iqb
        !> The tolerance factor is used to reduce the size of the Vc-diagonal product basis
        real(dp),     intent(in)                :: barcevtol
        !> The format in which the files format has been saved
        character(len=*), intent(in)            :: gw_io_format  
        !> The original matrix (the third index is the one of the matrices set)
        complex(dp), intent(in)                 :: matrix_qa(:,:,:)
        !> The new matrix (the third index is the one of the matrices set)
        complex(dp), allocatable, intent(out)   :: matrix_qb(:,:,:)

        ! Rotation matrices in MB and bare Coulomb eigenbasis
        complex(dp), allocatable :: rotmat_mb(:,:), rotmat_cb(:,:)

        ! Temporary for matrix-matrix product
        complex(dp), allocatable :: MatrixRot(:,:)

        ! Timings
        real(dp) :: time_start, time_end

        ! Index
        integer(i32) :: imat

        ! Start the timer
        call timesec(time_start)

        ! Compute the rotation in the mixed basis
        call compute_rotation_matrix_mb(iqa, iqb, kset, kqset, Gqset, Gset, gw_io_format, rotmat_mb)

        ! Basis change the rotation matrix to the Coulomb eigenbasis
        call rotmatmb_2_bare_coulomb_eigenbasis(iqa, iqb, kqset, barcevtol, gw_io_format, rotmat_mb, rotmat_cb)

        ! Free the rotation in mixed basis as it is not needed anymore
        deallocate(rotmat_mb)

        ! Rotate the original matrix from qa to qb
        allocate(MatrixRot, mold=matrix_qa(:,:,1))
        allocate(matrix_qb, mold=matrix_qa)

        ! $M(\mathbf{q}_b) = R^{\dagger}\cdot M(\mathbf{q}_a) \cdot R$
        ! Notice all are in their Coulomb eigenbasis
        do imat = 1, size(matrix_qa, 3)
            call matrix_multiply(matrix_qa(:,:,imat), rotmat_cb, MatrixRot)
            call matrix_multiply(rotmat_cb, MatrixRot, matrix_qb(:,:,imat), 'c', 'n')
        end do
        
        ! Add Timing
        call timesec(time_end)
        time_rotmb = time_rotmb + (time_end - time_start)

    end subroutine rotate_from_qa_to_qb_matrices

end module mod_gw_symmetry
