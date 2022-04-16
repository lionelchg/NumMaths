module mod_ode_analytic_brusselator

use, intrinsic :: iso_c_binding
use mod_prec

implicit none

integer(c_long), parameter :: neq = 3

integer :: n_problems

real(pr), parameter :: a  = 1.2d0
real(pr), parameter :: b  = 2.5d0
real(pr), parameter :: ep = 1.0d-5

real(pr), dimension(:,:), allocatable :: Jbuf_bruss



contains

integer(c_int) function RhsFn_bruss(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='RhsFn_bruss')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(pr), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr),    value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(pr), pointer :: yvec_1d(:)
    real(pr), pointer :: fvec_1d(:)
    real(pr), pointer :: yvec(:, :)
    real(pr), pointer :: fvec(:, :)
    integer :: i

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec_1d => FN_VGetArrayPointer(sunvec_y)
    fvec_1d => FN_VGetArrayPointer(sunvec_f)

    ! pointer bounds remapping
    yvec(1:neq, 1:n_problems) => yvec_1d
    fvec(1:neq, 1:n_problems) => fvec_1d

    ! fill RHS vector
    do i = 1, n_problems
        fvec(1, i) = a  -  (yvec(3, i) + 1.0d0) * yvec(1, i)  +  yvec(2, i) * yvec(1, i) * yvec(1, i)
        fvec(2, i) = yvec(3, i) * yvec(1, i)  -  yvec(2, i) * yvec(1, i) * yvec(1, i)
        fvec(3, i) = (b-yvec(3, i))/ep - yvec(3, i) * yvec(1, i)
    end do
    ! return success
    ierr = 0
    return

end function RhsFn_bruss

integer(c_int) function JacFn_bruss(tn, sunvec_y, sunvec_f, sunmat_J, &
       user_data, tmp1, tmp2, tmp3) &
       result(ierr) bind(C,name='JacFn_bruss')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsunmatrix_dense_mod
    use fsundials_matrix_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(pr), value :: tn                     ! current time
    type(N_Vector)        :: sunvec_y         ! current solution N_Vector
    type(N_Vector)        :: sunvec_f         ! current rhs N_Vector
    type(SUNMatrix)       :: sunmat_J         ! Jacobian SUNMatrix
    type(c_ptr), value    :: user_data        ! user-defined data
    type(N_Vector)        :: tmp1, tmp2, tmp3 ! workspace N_Vectors

    ! pointers to data in SUNDIALS vectors
    real(pr), pointer :: yvec_1d(:)
    real(pr), pointer :: yvec(:, :)
    integer(c_long) :: loc, np, icol, istart

    ! pointer to data in SUNDIALS matrix
    real(pr), pointer :: Jmat_vec(:)
    real(pr), pointer :: Jmat(:,:)
    real(pr), pointer :: Jcol(:)

    !======= Internals ============

    ! get data array from SUNDIALS vector
    yvec_1d => FN_VGetArrayPointer(sunvec_y)
    yvec(1:neq, 1:n_problems) => yvec_1d

    ! get data arrays from SUNDIALS vectors
    Jmat_vec => FSUNDenseMatrix_Data(sunmat_J)

    ! Set the jacobian to zero
    Jmat_vec(:) = 0.0D0

    ! Jacobian as a matrix for fortran
    Jmat(0:neq*n_problems-1, 0:neq*n_problems-1) => Jmat_vec

    ! Iterate over the columns of the Jacobian
    do np = 1, n_problems
        ! We need a pointer to the first column of the nonzeroblock of the jacobian for the current problem np
        icol = (np - 1) * neq
        ! Pointer to the ith column of the Jacobian matrix. WARNING: 0-indexing for the pointer
        Jcol(0:neq*n_problems-1) => FSUNDenseMatrix_Column(sunmat_J, icol)
        ! Now we set the correct subset of the column, which starts at
        istart = (np - 1) * neq
        Jcol(istart:istart+neq-1) = [ -(yvec(3, np)+1.0d0) + 2.0d0*yvec(1, np)*yvec(2, np), yvec(3, np) - 2.0d0*yvec(1, np)*yvec(2, np), -yvec(3, np)]

        ! Second column of the block
        Jcol(0:neq*n_problems-1) => FSUNDenseMatrix_Column(sunmat_J, icol+1)
        Jcol(istart:istart+neq-1) = [ yvec(1, np)*yvec(1, np), -yvec(1, np)*yvec(1, np), 0.0d0]

        ! Third column
        Jcol(0:neq*n_problems-1) => FSUNDenseMatrix_Column(sunmat_J, icol+2)
        Jcol(istart:istart+neq-1) = [ -yvec(1, np), yvec(1, np), -1.0d0/ep - yvec(1, np)]
    end do

    ! do np = 0, n_problems*neq - 1
    !     write (*, '(100(ES12.5,1X))') Jmat(:, np)
    ! end do
    ! call exit(1)

    ! return success
    ierr = 0
    return

  end function JacFn_bruss



end module