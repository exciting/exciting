!> Submodule implementing procedures from mod_APW_LO
submodule(mod_apw_lo) smod_apw_lo

contains
        !> Save APW and local-orbital data 
        module subroutine save_apwlo()

            use mod_eigensystem, only: oalo, ololo

            implicit none

            integer :: un

            open(newunit=un, file="apwlo.basis", form='unformatted', status='replace', action='write')

            
            !> APW data                                                        
            
            write(un) apword
            write(un) apwordmax
            write(un) apwe0
            write(un) apwdm
            write(un) apwn
            write(un) apwve

            !> apwe
            write(un) logical(allocated(apwe),kind=i32)
            if (allocated(apwe)) then
                write(un) lbound(apwe,kind=i32)
                write(un) ubound(apwe,kind=i32)
                write(un) apwe
            end if

            !> apwfr
            write(un) logical(allocated(apwfr),kind=i32)
            if (allocated(apwfr)) then
                write(un) lbound(apwfr,kind=i32)
                write(un) ubound(apwfr,kind=i32)
                write(un) apwfr
            end if

            !> apwdfr
            write(un) logical(allocated(apwdfr),kind=i32)
                if (allocated(apwdfr)) then
                write(un) lbound(apwdfr,kind=i32)
                write(un) ubound(apwdfr,kind=i32)
                write(un) apwdfr
            end if

            
            !> Local-orbital data

            write(un) nlorb
            write(un) nlomax
            write(un) nlotot
            write(un) lorbord
            write(un) lorbl
            write(un) lorbn
            write(un) lorbk
            write(un) wfkappa
            write(un) lolmax
            write(un) lolmmax
            write(un) lorbe0
            write(un) lorbdm
            write(un) lorbve
            write(un) lorbwfproj

            !> lorbe
            write(un) logical(allocated(lorbe),kind=i32)
            if (allocated(lorbe)) then
                write(un) lbound(lorbe,kind=i32)
                write(un) ubound(lorbe,kind=i32)
                write(un) lorbe
            end if

            !> lofr
            write(un) logical(allocated(lofr),kind=i32)
            if (allocated(lofr)) then
                write(un) lbound(lofr,kind=i32)
                write(un) ubound(lofr,kind=i32)
                write(un) lofr
            end if

            
            ! Miscellaneous                 
            write(un) mine0

            ! Overlap integrals
            write(un) logical(allocated(oalo),kind=i32)
            if (allocated(oalo)) then
                write(un) lbound(oalo,kind=i32)
                write(un) ubound(oalo,kind=i32)
                write(un) oalo
            end if

            write(un) logical(allocated(ololo),kind=i32)
            if (allocated(ololo)) then
                write(un) lbound(ololo,kind=i32)
                write(un) ubound(ololo,kind=i32)
                write(un) ololo
            end if

            close(un)

        end subroutine save_apwlo

        !> Load APW and local-orbital data 
        module subroutine load_apwlo()

            use mod_eigensystem, only: oalo, ololo

            implicit none

            integer :: un

            logical(i32) :: is_allocated

            integer(i32) :: lb3(3), ub3(3)
            integer(i32) :: lb4(4), ub4(4)
            integer(i32) :: lb5(5), ub5(5)

            open(newunit=un, file="apwlo.basis", form='unformatted', status='old', action='read')
            
            !> APW data                                                        !

            read(un) apword
            read(un) apwordmax
            read(un) apwe0
            read(un) apwdm
            read(un) apwn
            read(un) apwve

            !> apwe
            read(un) is_allocated
            if (is_allocated) then

                read(un) lb3
                read(un) ub3

                if (allocated(apwe)) deallocate(apwe)

                allocate(apwe( &
                        lb3(1):ub3(1), &
                        lb3(2):ub3(2), &
                        lb3(3):ub3(3)))

                read(un) apwe

            end if

            !> apwfr
            read(un) is_allocated
            if (is_allocated) then

                read(un) lb5
                read(un) ub5

                if (allocated(apwfr)) deallocate(apwfr)

                allocate(apwfr( &
                        lb5(1):ub5(1), &
                        lb5(2):ub5(2), &
                        lb5(3):ub5(3), &
                        lb5(4):ub5(4), &
                        lb5(5):ub5(5)))

                read(un) apwfr

            end if

            !> apwdfr
            read(un) is_allocated
            if (is_allocated) then

                read(un) lb3
                read(un) ub3

                if (allocated(apwdfr)) deallocate(apwdfr)

                allocate(apwdfr( &
                        lb3(1):ub3(1), &
                        lb3(2):ub3(2), &
                        lb3(3):ub3(3)))

                read(un) apwdfr

            end if

            
            !> Local-orbital data
            

            read(un) nlorb
            read(un) nlomax
            read(un) nlotot
            read(un) lorbord
            read(un) lorbl
            read(un) lorbn
            read(un) lorbk
            read(un) wfkappa
            read(un) lolmax
            read(un) lolmmax
            read(un) lorbe0
            read(un) lorbdm
            read(un) lorbve
            read(un) lorbwfproj

            !> lorbe
            read(un) is_allocated
            if (is_allocated) then

                read(un) lb3
                read(un) ub3

                if (allocated(lorbe)) deallocate(lorbe)

                allocate(lorbe( &
                        lb3(1):ub3(1), &
                        lb3(2):ub3(2), &
                        lb3(3):ub3(3)))

                read(un) lorbe

            end if

            !> lofr
            read(un) is_allocated
            if (is_allocated) then

                read(un) lb4
                read(un) ub4

                if (allocated(lofr)) deallocate(lofr)

                allocate(lofr( &
                        lb4(1):ub4(1), &
                        lb4(2):ub4(2), &
                        lb4(3):ub4(3), &
                        lb4(4):ub4(4)))

                read(un) lofr

            end if

            
            !> Miscellaneous
            read(un) mine0

            ! Overlap integrals
            read(un) is_allocated
            if (is_allocated) then
                
                read(un) lb3
                read(un) ub3
                
                if (allocated(oalo)) deallocate(oalo)

                allocate(oalo( &
                        lb3(1):ub3(1), &
                        lb3(2):ub3(2), &
                        lb3(3):ub3(3)))

                read(un) oalo

            end if

            read(un) is_allocated
            if (is_allocated) then
                
                read(un) lb3
                read(un) ub3
                
                if (allocated(ololo)) deallocate(ololo)

                allocate(ololo( &
                        lb3(1):ub3(1), &
                        lb3(2):ub3(2), &
                        lb3(3):ub3(3)))

                read(un) ololo

            end if

            close(un)

        end subroutine load_apwlo

end submodule smod_apw_lo