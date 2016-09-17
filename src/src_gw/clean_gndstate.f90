
subroutine clean_gndstate

    use modmain
    
    !if (allocated()) deallocate()
    
    !----------------
    ! init0
    !----------------
   
    if (allocated(rhomt)) deallocate(rhomt)
    if (allocated(rhoir)) deallocate(rhoir)
    if (allocated(magmt)) deallocate(magmt)
    if (allocated(magir)) deallocate(magir)
    if (allocated(vclmt)) deallocate(vclmt)
    if (allocated(vclir)) deallocate(vclir)
    if (allocated(vxcmt)) deallocate(vxcmt)
    if (allocated(vxcir)) deallocate(vxcir)      
    if (allocated(bxcmt)) deallocate(bxcmt)
    if (allocated(bxcir)) deallocate(bxcir)
    if (allocated(exmt)) deallocate(exmt)
    if (allocated(exir)) deallocate(exir)
    if (allocated(ecmt)) deallocate(ecmt)
    if (allocated(ecir)) deallocate(ecir)
    if (allocated(veffmt)) deallocate(veffmt)
    if (allocated(veffir)) deallocate(veffir)
    if (allocated(veffig)) deallocate(veffig)
    if (allocated(chgmt)) deallocate(chgmt)
    if (allocated(mommt)) deallocate(mommt)
    
    if (allocated(forcehf)) deallocate(forcehf)
    if (allocated(forcecr)) deallocate(forcecr)
    if (allocated(forceibs)) deallocate(forceibs)

    !----------------
    ! init1
    !----------------
    if (allocated(apwe)) deallocate(apwe)
    if (allocated(lorbe)) deallocate(lorbe)
    
    if (allocated(oalo)) deallocate(oalo)
    if (allocated(ololo)) deallocate(ololo)
    if (allocated(haa)) deallocate(haa)
    if (allocated(hloa)) deallocate(hloa)
    if (allocated(hlolo)) deallocate (hlolo)
    if (allocated(gntyry)) deallocate (gntyry)
    
    
end subroutine
