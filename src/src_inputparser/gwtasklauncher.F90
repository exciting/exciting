
subroutine gwtasklauncher()
    use modinput
    use inputdom
    use mod_selfconsistent_gw, only: is_gw_selfconsistent_flavour, qsgw

    if (.not. is_gw_selfconsistent_flavour(qsgw)) call rereadinput()
    call gw_main()

    return
end subroutine gwtasklauncher
