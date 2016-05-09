!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine idfgather
    Use modmain
    Use modinput
    Use modxs
#ifdef TETRA    
    Use modtetra
#endif    
    Use modmpi
    Use m_filedel
    Use m_getunit
    Use m_genfilname
    Implicit None
    ! local variables
    Character (*), Parameter :: thisnam = 'idfgather'
    Character (256) :: filnam, filnam2
    Integer :: n, m, iq, iw, iproc, recl, nc, oct1, oct2, octl, octu
    Logical :: tq0
    Complex (8), Allocatable :: mdf1 (:)
    Logical, External :: tqgamma
    Allocate (mdf1(nwdf))
    Inquire (IoLength=Recl) mdf1 (1)
    Call getunit (unit1)
    ! loop over q-points
    Do iq = 1, nqpt
        tq0 = tqgamma (iq)
        ! number of components (3 for q=0)
        nc = 1
        If (tq0) nc = 3
        ! matrix size for local field effects
        n = ngq (iq)
        ! calculate k+q and G+k+q related variables
        Call init1offs (qvkloff(1, iq))
        Do m = 1, n, Max (n-1, 1)
            ! loop over longitudinal components for optics
            Do oct1 = 1, nc
                If (input%xs%dfoffdiag) Then
                    octl = 1
                    octu = nc
                Else
                    octl = oct1
                    octu = oct1
                End If
                Do oct2 = octl, octu
                    Do iproc = 0, procs - 1
                        wpari = firstofset (iproc, nwdf)
                        wparf = lastofset (iproc, nwdf)
                        ! filename for proc
                        Call genfilname (basename='IDF', bzsampl=bzsampl, &
                        acont=input%xs%tddft%acont, nar= .Not. &
                        input%xs%tddft%aresdf, nlf=(m == 1), &
                        fxctype=input%xs%tddft%fxctypenumber, tq0=tq0, &
                        oc1=oct1, oc2=oct2, iqmt=iq, procs=procs, &
                        rank=iproc, filnam=filnam2)
                        if(rank.eq.iproc)then
                            Open (unit1, File=trim(filnam2), Form='unformatted'&
                            , Action='read', Status='old', Access='direct', &
                            Recl=Recl)
                            Do iw = wpari, wparf
                                Read (unit1, Rec=iw-wpari+1) mdf1 (iw)
                            End Do
                            Close (unit1)
                        endif
#ifdef MPI
                        call mpi_bcast(mdf1(wpari),wparf-wpari+1,MPI_DOUBLE_COMPLEX,iproc,MPI_COMM_WORLD,ierr)
#endif
                    End Do ! iproc
                    ! write to file
                    Call genfilname (basename='IDF', bzsampl=bzsampl, &
                    acont=input%xs%tddft%acont, nar= .Not. &
                    input%xs%tddft%aresdf, nlf=(m == 1), &
                    fxctype=input%xs%tddft%fxctypenumber, tq0=tq0, &
                    oc1=oct1, oc2=oct2, iqmt=iq, filnam=filnam)
                    if(rank.eq.0.or.(firstinnode.and. .not.input%sharedfs))then
                        Open (unit1, File=trim(filnam), Form='unformatted', &
                        Action='write', Status='replace', Access='direct', &
                        Recl=Recl)
                        Do iw = 1, nwdf
                            Write (unit1, Rec=iw) mdf1 (iw)
                        End Do
                        Close (unit1)
                    ENDIF
                    ! remove partial files
                    Do iproc = 0, procs - 1
                        if(iproc.eq.rank)then
                            Call genfilname (basename='IDF', bzsampl=bzsampl, &
                            acont=input%xs%tddft%acont, nar= .Not. &
                            input%xs%tddft%aresdf, nlf=(m .Eq. 1), &
                            fxctype=input%xs%tddft%fxctypenumber, tq0=tq0, &
                            oc1=oct1, oc2=oct2, iqmt=iq, procs=procs, &
                            rank=iproc, filnam=filnam2)
                            Call filedel (trim(filnam2))
                        endif
                    End Do
                ! end loop over optical components
                End Do
            End Do
        End Do ! m
        Write (unitout, '(a, i8)') 'Info(' // thisnam // '): inverse d&
        &ielectric function gathered for q - point:'        , iq
    End Do
    Deallocate (mdf1)
End Subroutine idfgather
