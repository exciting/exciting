!
!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!

Subroutine idfq(iq)
    Use modinput, Only: input
    Use modxs, Only: ivgigq, wpari, wparf, ngq, ikmapikq, istocc0, istunocc0, &
                    & isto0, isto, istu0, istu, ivgmt, unitout, bzsampl, unit1, &
                    & istocc, istunocc, nwdf, tfxcbse, qvkloff
    Use modfxcifc, Only: fxcifc
    Use constants, Only: zzero
    Use mod_misc, Only: task
#ifdef TETRA
    Use modtetra
#endif
    Use modmpi, Only: procs, rank, terminate
    Use m_genwgrid, Only: genwgrid
    Use m_dyson, Only: dyson
    Use m_dysonsym, Only: dysonsym
    Use m_getx0, Only: getx0
    Use m_genfilname, Only: genfilname
    Implicit None

    ! Arguments
    Integer, Intent(In) :: iq

    ! Local variables
    Character(*), Parameter :: thisnam = 'idfq'
    Character(256) :: filnam, filnam2
    Complex(8), Allocatable :: chi0(:, :), fxc(:, :), idf(:, :), mdf1(:), w(:)
    Complex(8), Allocatable :: chi0hd(:), chi0wg(:, :, :), chi0h(:, :)
    Integer :: n, m, recl, j, iw, wi, wf, nc, oct1, oct2, octl, octu, igmt
    Logical :: tq0
    Integer, External :: l2int
    Logical, External :: tqgamma

    ! Sampling type for Brillouin zone sampling
    bzsampl = l2int(input%xs%tetra%tetradf)
    tq0 = tqgamma(iq)

    ! Number of components (3 for q=0)
    nc = 1
    If (tq0) Then
        nc = 3
    End If

    ! Limits for w-points
    wi = wpari
    wf = wparf

    ! Matrix size for local field effects
    n = ngq(iq)
    Allocate(chi0(n, n), fxc(n, n), idf(n, n), w(nwdf), mdf1(nwdf), chi0hd(nwdf))
    Allocate(chi0wg(n, 2, 3), chi0h(3, 3))
    fxc = zzero

    ! Filename for response function file
    Call genfilname(basename='X0', &
        asc=.False., bzsampl=bzsampl, acont=input%xs%tddft%acont, &
        nar=.Not.input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, &
        markfxcbse=tfxcbse, iqmt=iq, filnam=filnam)
    Call genfilname(iqmt=iq, setfilext=.True.)
    Call init1offs(qvkloff(1, iq))

    ! Find highest (partially) occupied and lowest (partially) unoccupied states
    Call findocclims(iq, ikmapikq(:, iq), istocc0, istunocc0, isto0, isto, istu0, istu)
    istunocc = istunocc0
    istocc = istocc0

    ! Find limits for band combinations
    Call ematbdcmbs(input%xs%emattype)

    ! Generate energy grid
    Call genwgrid(nwdf, input%xs%energywindow%intv, input%xs%tddft%acont, 0.d0, w_cmplx=w)

    ! Record length
    Inquire(IoLength=Recl) mdf1(1)

    igmt = ivgigq(ivgmt(1, iq), ivgmt(2, iq), ivgmt(3, iq), iq)
    If (igmt .Gt. n) Then
        Write(*,*)
        Write(*, '("Error(", a, "): G-vector index for momentum transfer out of range: ", i8)') trim(thisnam), igmt
        Write(*,*)
        Call terminate
    End If

    If ((igmt .Ne. 1) .And. (iw .Eq. wi)) Then
        Write(unitout,*)
        Write(unitout, '("Info(", a, "): non-zero G-vector Fourier component for momentum transfer:")') trim(thisnam)
        Write(unitout, '(" G-vector number         :", i8)') igmt
        Write(unitout, '(" G-vector (latt. coords.):", 3i8)') ivgmt(:, iq)
        Write(unitout,*)
    End If

    ! Neglect/include local field effects
    Do m = 1, n, Max(n-1, 1)
        Select Case(input%xs%tddft%fxctypenumber)
        Case(5)
            ! fxc kernel: ALDA
            ! The ALDA kernel does not depend on q in principle, but the G-mesh
            ! depends through its cutoff for G+q on q. It is independent of w.
            Call fxcifc(input%xs%tddft%fxctypenumber, iq=iq, ng=m, fxcg=fxc)
            ! Add symmetrized Coulomb potential (is equal to unity matrix)
            If (m .Eq. 1) Then
                fxc(igmt, igmt) = fxc(igmt, igmt) + 1.d0
            Else
                forall(j=1:m) fxc(j,j) = fxc(j,j) + 1.d0
            End If
        End Select

        ! Loop over longitudinal components for optics
        Do oct1 = 1, nc
            If (input%xs%dfoffdiag) Then
                octl = 1
                octu = nc
            Else
                octl = oct1
                octu = oct1
            End If

            Do oct2 = octl, octu
                ! Filename for output file
                Call genfilname(basename='IDF', &
                    asc=.False., bzsampl=bzsampl, acont=input%xs%tddft%acont, &
                    nar=.Not.input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m .Eq. 1), &
                    fxctype=input%xs%tddft%fxctypenumber, &
                    tq0=tq0, oc1=oct1, oc2=oct2, iqmt=iq, procs=procs, rank=rank, filnam=filnam2)
                Open(newunit=unit1, File=trim(filnam2), Form='unformatted', Action='write', Access='direct', Recl=Recl)

                Select Case(input%xs%tddft%fxctypenumber)
                Case(9:11)
                    ! fxc kernels: BO, BO_SCALAR, RBO
                    Call getx0(tq0, iq, 1, trim(filnam), '', chi0, chi0wg, chi0h)
                    If (tq0) Then
                        chi0(1, 1) = chi0h(oct1, oct2)
                        If (m .Gt. 1) Then
                            chi0(1, 2:) = chi0wg(2:, 1, oct1)
                            chi0(2:, 1) = chi0wg(2:, 2, oct2)
                        End If
                    End If
                    Call fxcifc(input%xs%tddft%fxctypenumber, ng=m, chim=chi0, fxcg=fxc)
                    ! Add symmetrized Coulomb potential (is equal to unity matrix)
                    If (m .Eq. 1) Then
                        fxc(igmt, igmt) = fxc(igmt, igmt) + 1.d0
                    Else
                        forall(j=1:m) fxc(j,j) = fxc(j,j) + 1.d0
                    End If
                End Select

                Do iw = wi, wf
                    Call chkpt(6, (/ task, iq, m, oct1, oct2, iw /), 'task, q - point index, loc. field., opt. comp. 1, opt. comp. 2, w - point; Dyson equation')

                    ! Read Kohn-Sham response function
                    Call getx0(tq0, iq, iw, trim(filnam), '', chi0, chi0wg, chi0h)

                    ! Assign components to main matrix for q=0
                    If (tq0) Then
                        ! Head
                        chi0(1, 1) = chi0h(oct1, oct2)
                        ! Wings
                        If (m .Gt. 1) Then
                            chi0(1, 2:) = chi0wg(2:, 1, oct1)
                            chi0(2:, 1) = chi0wg(2:, 2, oct2)
                        End If
                    End If

                    ! Generate xc-kernel
                    Select Case(input%xs%tddft%fxctypenumber)
                    Case(0, 1, 2, 3, 4)
                        ! fxc kernels: RPA, LRCstatic_NLF, LRCstatic, LRCdyn_NLF, LRCdyn
                        Call fxcifc(input%xs%tddft%fxctypenumber, &
                            ng=m, iw=iw, w=w(iw), &
                            alrc=input%xs%tddft%alphalrc, &
                            alrcd=input%xs%tddft%alphalrcdyn, &
                            blrcd=input%xs%tddft%betalrcdyn, &
                            fxcg=fxc)
                    Case(7, 8)
                        ! fxc kernels: MB1_NLF, MB1
                        Call fxcifc(input%xs%tddft%fxctypenumber, &
                            oct=oct1, ng=m, iw=iw, w=w(iw), &
                            alrc=input%xs%tddft%alphalrc, &
                            alrcd=input%xs%tddft%alphalrcdyn, &
                            blrcd=input%xs%tddft%betalrcdyn, &
                            fxcg=fxc)
                    End Select                   

                    ! Solve Dyson's equation for the interacting response function
                    Select Case(input%xs%tddft%fxctypenumber)
                    Case(0, 1, 2, 3, 4)
                        ! fxc kernels: RPA, LRCstatic_NLF, LRCstatic, LRCdyn_NLF, LRCdyn
                        ! Add symmetrized Coulomb potential (is equal to unity matrix)
                        If (m .Eq. 1) Then
                            fxc(igmt, igmt) = fxc(igmt, igmt) + 1.d0
                        Else
                            forall(j=1:m) fxc(j,j) = fxc(j,j) + 1.d0
                        End If
                        Call dyson(n, chi0, fxc, idf)
                    Case(5, 9:11)
                        ! fxc kernels: ALDA, BO, BO_SCALAR, RBO
                        Call dyson(n, chi0, fxc, idf)
                    Case(7, 8)
                        ! fxc kernels: MB1_NLF, MB1
                        ! We do not expect the kernel to contain the symmetrized
                        ! Coulomb potential here, the kernel here is expected to be
                        ! multiplied with the KS response function from both sides.
                        ! [F. Sottile, PRL 2003]
                        Call dysonsym(n, chi0, fxc, idf)
                    End Select

                    ! Symmetrized inverse dielectric function (add one)
                    If (m .Eq. 1) Then
                        idf(igmt, igmt) = idf(igmt, igmt) + 1.d0
                    Else
                        forall(j=1:m) idf(j,j) = idf(j,j) + 1.d0
                    End If

                    ! Adler-Wiser treatment of macroscopic dielectric function
                    mdf1(iw) = 1.d0 / idf(igmt, igmt)

                    ! Mimic zero Kronecker delta in case of off-diagonal tensor components (?)
                    If ((m .Eq. 1) .And. (oct1 .Ne. oct2)) mdf1(iw) = mdf1(iw) - 1.d0

                    ! Write macroscopic dielectric function to file
                    Write(unit1, Rec=iw-wi+1) mdf1(iw)
                End Do ! iw
                Close(unit1)
            End Do ! oct2
        End Do ! oct1
    End Do ! m

     ! Deallocate
    Deallocate(chi0, chi0wg, chi0h, fxc, idf, mdf1, w, chi0hd)
End Subroutine idfq
