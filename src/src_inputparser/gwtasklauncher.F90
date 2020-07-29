
subroutine gwtasklauncher()
    use modinput
    use inputdom

    call rereadinput()
    call gw_main()

    return
end subroutine gwtasklauncher
