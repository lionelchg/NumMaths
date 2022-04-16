module mod_implicit_solver
    use, intrinsic :: iso_c_binding
    use mod_prec

    use fcvode_mod                      ! Fortran interface to CVODE
    use fnvector_serial_mod             ! Fortran interface to serial N_Vector
    use fsunmatrix_dense_mod            ! Fortran interface to dense SUNMatrix
    use fsunlinsol_dense_mod            ! Fortran interface to dense SUNLinearSolver
    use fsundials_linearsolver_mod      ! Fortran interface to generic SUNLinearSolver
    use fsundials_matrix_mod            ! Fortran interface to generic SUNMatrix
    use fsundials_nvector_mod           ! Fortran interface to generic N_Vector

    implicit none


    real(pr) :: rtol, atol           ! Tolerances
    integer(c_long) :: neq

    type(c_ptr) :: cvode_mem        ! CVODE memory

    real(pr), dimension(:), allocatable :: ydata
    type(N_Vector),        pointer :: ysun      ! Sundials vector
    type(SUNMatrix),       pointer :: Asun      ! Sundials matrix
    type(SUNLinearSolver), pointer :: LSsun     ! Sundials linear solver


    contains


    subroutine init ( neq_in, rtol_in, atol_in, system )
        use mod_ode_analytic_boeuf, only: RhsFn
        use mod_ode_analytic_brusselator, only: RhsFn_bruss, JacFn_bruss, Jbuf_bruss, n_problems

        ! IN/OUT
        integer(c_long) :: neq_in
        real(pr) :: rtol_in, atol_in
        integer :: system

        ! LOCAL
        integer :: ierr
        real(pr) :: dummy_real

        ! Set attributes
        rtol = rtol_in
        atol = atol_in
        neq = neq_in
        dummy_real = 0.0D0

        write (*, '(A, I0)') "allocating for neq = ", neq

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
        if ( system == 0 ) then
            ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), dummy_real, ysun)
        else if ( system == 1 .or. system == 2) then
            ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn_bruss), dummy_real, ysun)
        end if
        if (ierr /= 0) then
            print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! More internal steps
        ierr = FCVodeSetMaxNumSteps(cvode_mem, 1000 )

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
        if ( system == 2 ) then
            ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn_bruss))
            ! Allocate jac
        end if

        ! Set pointer to user_data array
        ! n_user_data = user_data_size
        ! allocate ( user_data(n_user_data) )
        ! user_data = 0.0D0
        ! ierr = FCVodeSetUserData(cvode_mem, c_loc(user_data))

    end subroutine

    subroutine solve ( y, tout )
        real(c_double), dimension(neq), intent(inout) :: y
        real(c_double), intent(in) :: tout
        real(c_double):: tcur(1)

        integer :: ierr
        real(8), pointer :: y_loc(:)
        integer(c_long), dimension(1) :: nsteps, nrhsevals, njacevals

        ydata(:) = y(:)

        ierr = FCVodeReInit(cvode_mem, 0.0D0, ysun)
        if (ierr /= 0) then
            print *, 'Error in FCVodeReInit, ierr = ', ierr, '; halting'
            stop 1
        end if

        ierr = FCVode(cvode_mem, tout, ysun, tcur(1), CV_NORMAL)
        if (ierr /= 0) then
            print *, 'Error in FCVODE, ierr = ', ierr, '; halting'
            stop 1
        end if

        ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
        ierr = FCVodeGetNumRhsEvals(cvode_mem, nrhsevals)
        ierr = FCVodeGetNumJacEvals(cvode_mem, njacevals)

        ! write (*, '(3(A,I0,1X))') "num_steps=", nsteps, "num_rhs_evals=", nrhsevals, "num_jac_evals=", njacevals

        y(:) = ydata(:)

    end subroutine solve

    subroutine destroy ( )

        integer :: ierr

        call FCVodeFree(cvode_mem)
        call FN_VDestroy(ysun)

        deallocate ( ydata )

    end subroutine destroy

end module mod_implicit_solver