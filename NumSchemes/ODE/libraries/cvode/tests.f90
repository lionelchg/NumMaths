program main
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod

    implicit none

    integer, parameter :: N = 10

    type(N_Vector), pointer :: xnv, ynv, znv
    real(8)        :: x(N), y(N), z(N)

    ! Initialise x
    x(:) = 10.0d0
    y(:) = 5.0D0
    z(:) = 0.0D0

    ! Create nvectors
    xnv => FN_VMake_Serial(N, x)
    ynv => FN_VMake_Serial(N, y)
    znv => FN_VMake_Serial(N, z)

    call FN_VLinearSum(1.0d0, xnv, 1.0d0, ynv, znv)

    print *, z

    y(:) = 2.0D0

    call FN_VLinearSum(1.0d0, xnv, 1.0d0, ynv, znv)

    print *, z



end program main