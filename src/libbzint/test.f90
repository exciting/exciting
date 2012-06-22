
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

program test

integer(4):: i,j,ind

      do i=1,4
       do j=1,4
        ind=mod(j+i-2,4)+1
        write(*,*)i,j,ind
       enddo
      enddo
end
