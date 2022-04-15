module mod_ode_analytic
use, intrinsic :: iso_c_binding

implicit none

integer(c_long), parameter :: neq = 3

! ODE parameters
! real(8), parameter :: a  = 1.2d0
! real(8), parameter :: b  = 2.5d0
! real(8), parameter :: ep = 1.0d-5
! real(8) :: a = 1.2d0
! real(8) :: b = 2.5d0
! real(8) :: ep = 1.0d-5

real(8), pointer :: a_ptr
real(8), pointer :: b_ptr
real(8), pointer :: ep_ptr


contains

integer(c_int) function RhsFn(t, sunvec_y, sunvec_f, user_data) &
        result(ierr) bind(C, name="RhsFn")
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    ! IN/OUT
    real(c_double), value :: t
    type(N_Vector) :: sunvec_y
    type(N_Vector) :: sunvec_f
    type(N_Vector) :: test
    type(c_ptr), value :: user_data

    ! Pointers to data in Sundials vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)

    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)

    ! write (*, '(3(A, ES12.5))') "a = ", a_ptr, " b = ", b_ptr, " ep = ", ep_ptr

    ! Fill RHS vector
    fvec(1) = a_ptr  -  (yvec(3) + 1.0d0) * yvec(1)  +  yvec(2) * yvec(1) * yvec(1)
    fvec(2) = yvec(3) * yvec(1)  -  yvec(2) * yvec(1) * yvec(1)
    fvec(3) = (b_ptr-yvec(3))/ep_ptr - yvec(3) * yvec(1)

    ! Return success
    ierr = 0
    return

end function RhsFn

integer(c_int) function JacFn(t, sunvec_y, sunvec_f, sunmat_J, user_data, tmp1, tmp2, tmp3) &
        result(ierr) bind(C, name='JacFn')
    use fsundials_nvector_mod
    use fsunmatrix_dense_mod
    use fsundials_matrix_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: sunvec_y, sunvec_f
    type(SUNMatrix) :: sunmat_J
    type(c_ptr), value :: user_data
    type(N_Vector) :: tmp1, tmp2, tmp3

    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: Jmat(:)

    yvec => FN_VGetArrayPointer(sunvec_y)
    Jmat => FSUNDenseMatrix_Data(sunmat_J)

    Jmat = [-(yvec(3)+1.0d0) + 2.0d0*yvec(1)*yvec(2),  &
            yvec(3) - 2.0d0*yvec(1)*yvec(2), -yvec(3), &
            yvec(1)*yvec(1), -yvec(1)*yvec(1), 0.0d0,  &
            -yvec(1), yvec(1), -1.0d0/ep_ptr - yvec(1)]

    ierr = 0
    return
end function JacFn

end module mod_ode_analytic


! module mod_implicit_solver
! use, intrinsic :: iso_c_binding

! use fcvode_mod                      ! Fortran interface to CVODE
! use fnvector_serial_mod             ! Fortran interface to serial N_Vector
! use fsunmatrix_dense_mod            ! Fortran interface to dense SUNMatrix
! use fsunlinsol_dense_mod            ! Fortran interface to dense SUNLinearSolver
! use fsundials_linearsolver_mod      ! Fortran interface to generic SUNLinearSolver
! use fsundials_matrix_mod            ! Fortran interface to generic SUNMatrix
! use fsundials_nvector_mod           ! Fortran interface to generic N_Vector

! implicit none


! type implicit_solver_t
!     real(8) :: rtol, atol           ! Tolerances
!     integer(c_long) :: neq

!     type(c_ptr) :: cvode_mem        ! CVODE memory

!     real(8), dimension(:), allocatable :: ydata
!     type(N_Vector),        pointer :: ysun      ! Sundials vector
!     type(SUNMatrix),       pointer :: Asun      ! Sundials matrix
!     type(SUNLinearSolver), pointer :: LSsun     ! Sundials linear solver

! contains
!     procedure :: init
!     procedure :: solve
!     ! procedure :: info
!     procedure :: destroy
! end type implicit_solver_t


! contains


! subroutine init ( self, neq, rtol, atol )
!     use mod_ode_analytic, only: RhsFn

!     ! IN/OUT
!     class(implicit_solver_t) :: self
!     integer(c_long) :: neq
!     real(8) :: rtol, atol

!     ! LOCAL
!     integer :: ierr
!     real(8) :: dummy_real

!     print *, "allocating for neq = ", neq

!     ! Set attributes
!     self%rtol = rtol
!     self%atol = atol
!     self%neq = neq
!     dummy_real = 0.0D0

!     ! Create ydata Fortran vector
!     allocate ( self%ydata(neq) )

!     ! Create the y N_Vector around ydata
!     self%ysun => FN_VMake_Serial(neq, self%ydata)
!     if (.not. associated(self%ysun)) then
!         print *, 'ERROR: sunvec = NULL'
!         stop 1
!      end if

!     ! Create a dense matrix
!     self%Asun => FSUNDenseMatrix(neq, neq)
!     if (.not. associated(self%Asun)) then
!         print *, 'ERROR: sunmat = NULL'
!         stop 1
!      end if

!     ! Create a dense linear solver
!     self%LSsun => FSUNDenseLinearSolver(self%ysun, self%Asun)
!     if (.not. associated(self%LSsun)) then
!         print *, 'ERROR: sunlinsol = NULL'
!         stop 1
!      end if

!     ! Create CVODE memory with BDF method for stiff problems
!     self%cvode_mem = FCVodeCreate(CV_BDF)
!     if (.not. c_associated(self%cvode_mem)) then
!         print *, 'ERROR: cvode_mem = NULL'
!         stop 1
!      end if

!     ! Initialize CVode
!     ierr = FCVodeInit(self%cvode_mem, c_funloc(RhsFn), dummy_real, self%ysun)
!     if (ierr /= 0) then
!         print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
!         stop 1
!      end if

!     ! Set tolerances
!     ierr = FCVodeSStolerances(self%cvode_mem, self%rtol, self%atol)
!     if (ierr /= 0) then
!         print *, 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
!         stop 1
!      end if

!     ! Attach linear solver
!     ierr = FCVodeSetLinearSolver(self%cvode_mem, self%LSsun, self%Asun)
!     if (ierr /= 0) then
!         print *, 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
!         stop 1
!      end if

! end subroutine

! subroutine solve ( self, y, tout, tcur)
!     class(implicit_solver_t) :: self
!     real(c_double), dimension(self%neq), intent(inout) :: y
!     real(c_double), intent(in) :: tout
!     real(c_double):: tcur(1)

!     integer :: ierr
!     real(8), pointer :: ydata(:)

!     print *, "titi"
!     ! self%ydata(:) = y(:)
!     ydata => FN_VGetArrayPointer(self%ysun)
!     ydata(:) = y(:)

!     print *, "tutu"
!     ! ierr = FCVodeReInit(self%cvode_mem, tcur(1), self%ysun)
!     if (ierr /= 0) then
!         print *, 'Error in FCVodeReInit, ierr = ', ierr, '; halting'
!         stop 1
!      end if

!     print *, "toto"
!     ierr = FCVode(self%cvode_mem, tout, self%ysun, tcur(1), CV_NORMAL)
!     if (ierr /= 0) then
!         print *, 'Error in FCVODE, ierr = ', ierr, '; halting'
!         stop 1
!      end if

!     print *, "tata"
!     ! y(:) = self%ydata(:)
!     y = ydata

! end subroutine solve

! subroutine destroy ( self )
!     class(implicit_solver_t) :: self

!     integer :: ierr

!     call FCVodeFree(self%cvode_mem)
!     call FN_VDestroy(self%ysun)

! end subroutine destroy

! end module mod_implicit_solver

module mod_implicit_solver
    use, intrinsic :: iso_c_binding

    use fcvode_mod                      ! Fortran interface to CVODE
    use fnvector_serial_mod             ! Fortran interface to serial N_Vector
    use fsunmatrix_dense_mod            ! Fortran interface to dense SUNMatrix
    use fsunlinsol_dense_mod            ! Fortran interface to dense SUNLinearSolver
    use fsundials_linearsolver_mod      ! Fortran interface to generic SUNLinearSolver
    use fsundials_matrix_mod            ! Fortran interface to generic SUNMatrix
    use fsundials_nvector_mod           ! Fortran interface to generic N_Vector

    implicit none


    real(8) :: rtol, atol           ! Tolerances
    integer(c_long) :: neq

    type(c_ptr) :: cvode_mem        ! CVODE memory

    real(8), dimension(:), allocatable :: ydata
    type(N_Vector),        pointer :: ysun      ! Sundials vector
    type(SUNMatrix),       pointer :: Asun      ! Sundials matrix
    type(SUNLinearSolver), pointer :: LSsun     ! Sundials linear solver


    contains


    subroutine init ( neq_in, rtol_in, atol_in )
        use mod_ode_analytic, only: RhsFn, JacFn

        ! IN/OUT
        integer(c_long) :: neq_in
        real(8) :: rtol_in, atol_in

        ! LOCAL
        integer :: ierr
        real(8) :: dummy_real

        ! Set attributes
        rtol = rtol_in
        atol = atol_in
        neq = neq_in
        dummy_real = 0.0D0

        print *, "allocating for neq = ", neq

        ! Create ydata Fortran vector
        allocate ( ydata(neq) )
        ydata = 0.0D0

        ! Create the y N_Vector around ydata
        ysun => FN_VMake_Serial(neq, ydata)
        if (.not. associated(ysun)) then
            print *, 'ERROR: sunvec = NULL'
            stop 1
        end if

        ! Create a dense matrix
        Asun => FSUNDenseMatrix(neq, neq)
        if (.not. associated(Asun)) then
            print *, 'ERROR: sunmat = NULL'
            stop 1
        end if

        ! Create a dense linear solver
        LSsun => FSUNDenseLinearSolver(ysun, Asun)
        if (.not. associated(LSsun)) then
            print *, 'ERROR: sunlinsol = NULL'
            stop 1
        end if

        ! Create CVODE memory with BDF method for stiff problems
        cvode_mem = FCVodeCreate(CV_BDF)
        if (.not. c_associated(cvode_mem)) then
            print *, 'ERROR: cvode_mem = NULL'
            stop 1
        end if

        ! Initialize CVode
        ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), dummy_real, ysun)
        if (ierr /= 0) then
            print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! Set tolerances
        ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! Attach linear solver
        ierr = FCVodeSetLinearSolver(cvode_mem, LSsun, Asun)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! Set Jacobian routine
        ! ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))

        ! Set pointer to user_data array
        ! n_user_data = user_data_size
        ! allocate ( user_data(n_user_data) )
        ! user_data = 0.0D0
        ! ierr = FCVodeSetUserData(cvode_mem, c_loc(user_data))

    end subroutine

    subroutine solve ( y, tout, tcur )
        real(c_double), dimension(neq), intent(inout) :: y
        real(c_double), intent(in) :: tout
        real(c_double):: tcur(1)

        integer :: ierr
        real(8), pointer :: y_loc(:)

        ydata(:) = y(:)

        ierr = FCVodeReInit(cvode_mem, tcur(1), ysun)
        if (ierr /= 0) then
            print *, 'Error in FCVodeReInit, ierr = ', ierr, '; halting'
            stop 1
        end if

        ierr = FCVode(cvode_mem, tout, ysun, tcur(1), CV_NORMAL)
        if (ierr /= 0) then
            print *, 'Error in FCVODE, ierr = ', ierr, '; halting'
            stop 1
        end if

        y(:) = ydata(:)

    end subroutine solve

    subroutine destroy ( )

        integer :: ierr

        call FCVodeFree(cvode_mem)
        call FN_VDestroy(ysun)

        deallocate ( ydata )

    end subroutine destroy

end module mod_implicit_solver



program main
    use mod_implicit_solver, only: init, solve, destroy
    ! use mod_ode_analytic, only: neq, RhsFn, a, b, ep
    use mod_ode_analytic, only: neq, RhsFn, a_ptr, b_ptr, ep_ptr
    implicit none

    real(8) :: tstart, tend, rtol, atol, dtout, tout
    real(8) :: tcur(1)
    integer :: ierr, nout, outstep
    real(8), dimension(:), allocatable :: y
    real(8), target :: a, b, ep

    ! type(implicit_solver_t) :: solver

    allocate( y(neq) )

    tstart = 0.0D0
    tend = 10.0D0
    tcur = tstart
    tout = tstart
    dtout = (tend-tstart)/10.d0
    nout = ceiling(tend/dtout)
    rtol = 1.0D-5
    atol = 1.0D-10

    y(1) = 3.9D0
    y(2) = 1.1D0
    y(3) = 2.8D0

    a = 1.2d0
    b = 2.5d0
    ep = 1.0d-5

    a_ptr => a
    b_ptr => b
    ep_ptr => ep

    print *, "init"
    ! call solver%init(neq, rtol, atol)
    call init(neq, rtol, atol)

    ! start time stepping
    print *, '   '
    print *, 'Finished initialization, starting time steps'
    print *, '   '
    print *, '      t           u           v           w'
    print *, '----------------------------------------------------'
    print '(1x,4(es12.5,1x))', tcur, y(1), y(2), y(3)

    print *, "loop"
    do outstep = 1, nout
        tout = min(tout + dtout, tend)
        ! tout = tend
        ! call solver%solve(y, tout, tcur)
        call solve(y(1), tout, tcur)
        a = a + 1.0d-2
        b = b + 1.0d-2
        ep = ep + 1.0d-7

        print '(1x,4(es12.5,1x))', tcur, y(1), y(2), y(3)

    end do


    ! call solver%destroy()
    call destroy()

end program main
