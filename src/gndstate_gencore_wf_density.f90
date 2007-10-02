subroutine  gndstate_gencore_wf_density

    call gencore
     ! find the new linearisation energies
     call linengy
     ! write out the linearisation energies
     call writelinen
     ! generate the APW radial functions
     call genapwfr
     ! generate the local-orbital radial functions
     call genlofr
     ! compute the overlap radial integrals
     call olprad
     ! compute the Hamiltonian radial integrals
     call hmlrad
end subroutine	
