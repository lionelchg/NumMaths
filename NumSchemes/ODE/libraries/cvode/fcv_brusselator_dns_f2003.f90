! ------------------------------------------------------------------
! Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
!                Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! This test simulates a brusselator problem from chemical kinetics.
! This is an ODE system with 3 components, Y = [u,v,w], satisfying
! the equations,
!
!    du/dt = a - (w+1)*u + v*u^2
!    dv/dt = w*u - v*u^2
!    dw/dt = (b-w)/ep - w*u
!
! for t in the interval [0.0, 10.0], with initial conditions
! Y0 = [3.9, 1.1, 2.8], and parameter values a=1.2, b=2.5, and
! ep=1.0e-5
!
! Here, all three components exhibit a rapid transient change during
! the first 0.2 time units, followed by a slow and smooth evolution.
! ------------------------------------------------------------------

module ode_mod
    use, intrinsic :: iso_c_binding
    implicit none

    integer*8, parameter :: neq = 3

    real(8), parameter :: a = 1.2d0
    real(8), parameter :: b = 2.5d0
    real(8), parameter :: ep = 1.0d-5

contains

    ! ----------------------------------------------------------------
    ! RhsFn provides the right hand side function for the
    ! ODE: dy/dt = f(t,y)
    !
    ! Return values:
    !    0 = success,
    !    1 = recoverable error,
    !   -1 = non-recoverable error
    ! ----------------------------------------------------------------
    integer function RhsFn(tn, sunvec_y, sunvec_f, user_data) result(ierr) bind(C, name="RhsFn")
        use fsundials_nvector_mod

        real(8),      value :: tn ! current time
        type(N_Vector)      :: sunvec_y
        type(N_Vector)      :: sunvec_f
        type(c_ptr),  value :: user_data

        real(8), pointer :: yvec(:)
        real(8), pointer :: fvec(:)


        ! get data arrays from sundials vectors
        yvec => FN_VGetArrayPointer(sunvec_y)
        fvec => FN_VGetArrayPointer(sunvec_f)

        ! Fill RHS vector
        fvec(1) = a  -  (yvec(3) + 1.0d0) * yvec(1)  +  yvec(2) * yvec(1) * yvec(1)
        fvec(2) = yvec(3) * yvec(1)  -  yvec(2) * yvec(1) * yvec(1)
        fvec(3) = (b-yvec(3))/ep - yvec(3) * yvec(1)

        ! Return success
        ierr = 0
        return

    end function RhsFn

    ! ----------------------------------------------------------------
    ! JacFn: The Jacobian of the ODE hand side function J = df/dy
    !
    ! Return values:
    !    0 = success,
    !    1 = recoverable error,
    !   -1 = non-recoverable error
    ! ----------------------------------------------------------------
    integer function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, tmp1, tmp2, tmp3) &
    result(ierr) bind(C, name="JacFn")

        use fsundials_nvector_mod
        use fsunmatrix_dense_mod
        use fsundials_matrix_mod

        implicit none

        real(8), value      :: tn
        type(N_Vector)      :: sunvec_y
        type(N_Vector)      :: sunvec_f
        type(SUNMatrix)     :: sunmat_J
        type(c_ptr), value  :: user_data
        type(N_Vector)      :: tmp1, tmp2, tmp3  ! Workspace nvectors

        ! Pointer to data in sundials vectors
        real(8), pointer :: yvec(:)

        ! Pointer to data in sundials matrix
        real(8), pointer :: Jmat(:)


        ! Get data array from sundials vector
        yvec => FN_VGetArrayPointer(sunvec_y)

        ! Get data array from sundials matrix
        Jmat => FSUNDenseMatrix_Data(sunmat_J)

        ! Fill Jacobian Matrix
        Jmat = [                                         &
            - (yvec(3) + 1.0d0) + 2.0d0*yvec(1)*yvec(2), &
            yvec(3) - 2.0d0*yvec(1)*yvec(2), -yvec(3),   &
            yvec(1)*yvec(1), - yvec(1)*yvec(1), 0.0d0,   &
            - yvec(1), yvec(1), -1.0d0/ep - yvec(1)      &
        ]

        ! return success
        ierr = 0
        return

    end function JacFn

end module ode_mod


program main
    use, intrinsic :: iso_c_binding

    use fcvode_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsunlinsol_dense_mod
    use fsundials_linearsolver_mod
    use fsundials_matrix_mod
    use fsundials_nvector_mod
    use ode_mod

    implicit none

    real(8) :: tstart
    real(8) :: tend
    real(8) :: rtol, atol
    real(8) :: dtout
    real(8) :: tout
    real(8) :: tcur(1)

    integer :: ierr
    integer :: nout

    integer :: outstep

    type(c_ptr) :: cvode_mem
    type(N_Vector), pointer :: sunvec_y
    type(SUNMatrix), pointer :: sunmat_A
    type(SUNLinearSolver), pointer :: sunlinsol_LS

    real(8) :: yvec(neq)



    ! initialize ODE
    tstart = 0.0d0
    tend   = 10.0d0
    tcur   = tstart
    tout   = tstart
    dtout  = (tend-tstart)/10.d0
    nout   = ceiling(tend/dtout)

    ! initialize solution vector
    yvec(1) = 3.9d0
    yvec(2) = 1.1d0
    yvec(3) = 2.8d0

    ! Create sundials nvector
    sunvec_y => FN_VMake_Serial(neq, yvec)

    ! Create a dense matrix
    sunmat_A => FSUNDenseMatrix(neq, neq)

    ! Create a dense linear solver
    sunlinsol_LS => FSUNDenseLinearSolver(sunvec_y, sunmat_A)

    ! Create CVODE memory
    cvode_mem = FCVodeCreate(CV_BDF)

    ! Initialize CVODE
    ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)

    ! Set tolerances
    rtol = 1.0d-5
    atol = 1.0d-10

    ierr = FCVodeSStolerances(cvode_mem, rtol, atol)

    ! Attach linear solver
    ierr = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)

    ! Set Jacobian routine
    ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))

    ! start time stepping
    print *, '   '
    print *, 'Finished initialization, starting time steps'
    print *, '   '
    print *, '      t           u           v           w'
    print *, '----------------------------------------------------'
    print '(1x,4(es12.5,1x))', tcur, yvec(1), yvec(2), yvec(3)

    do outstep = 1, nout

        ! call CVODE
        tout = min(tout + dtout, tend)
        ierr = FCVode(cvode_mem, tout, sunvec_y, tcur, CV_NORMAL)

        ! output current solution
        print '(1x,4(es12.5,1x))', tcur, yvec(1), yvec(2), yvec(3)

    end do

    ! diagnostics output
    call CVodeStats(cvode_mem)

    ! clean up
    call FCVodeFree(cvode_mem)
    ierr = FSUNLinSolFree(sunlinsol_LS)
    call FSUNMatDestroy(sunmat_A)
    call FN_VDestroy(sunvec_y)


end program main

! ----------------------------------------------------------------
! CVodeStats
!
! Print CVODE statstics to standard out
! ----------------------------------------------------------------
subroutine CVodeStats(cvode_mem)

    !======= Inclusions ===========
    use iso_c_binding
    use fcvode_mod

    !======= Declarations =========
    implicit none

    type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

    integer(c_int)  :: ierr          ! error flag

    integer(c_long) :: nsteps(1)     ! num steps
    integer(c_long) :: nfevals(1)    ! num function evals
    integer(c_long) :: nlinsetups(1) ! num linear solver setups
    integer(c_long) :: netfails(1)   ! num error test fails

    integer(c_int)  :: qlast(1)      ! method order in last step
    integer(c_int)  :: qcur(1)       ! method order for next step

    real(c_double)  :: hinused(1)    ! initial step size
    real(c_double)  :: hlast(1)      ! last step size
    real(c_double)  :: hcur(1)       ! step size for next step
    real(c_double)  :: tcur(1)       ! internal time reached

    integer(c_long) :: nniters(1)    ! nonlinear solver iterations
    integer(c_long) :: nncfails(1)   ! nonlinear solver fails

    integer(c_long) :: njevals(1)    ! num Jacobian evaluations

    !======= Internals ============

    ! general solver statistics
    ierr = FCVodeGetIntegratorStats(cvode_mem, nsteps, nfevals, nlinsetups, &
         netfails, qlast, qcur, hinused, hlast, hcur, tcur)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetIntegratorStats, ierr = ', ierr, '; halting'
       stop 1
    end if

    ! nonlinear solver statistics
    ierr = FCVodeGetNonlinSolvStats(cvode_mem, nniters, nncfails)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNonlinSolvStats, ierr = ', ierr, '; halting'
       stop 1
    end if

    ! nonlinear solver statistics
    ierr = FCVodeGetNumJacEvals(cvode_mem, njevals)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumJacEvals, ierr = ', ierr, '; halting'
       stop 1
    end if

    print *, ' '
    print *, ' General Solver Stats:'
    print '(4x,A,i9)'    ,'Total internal steps taken =',nsteps
    print '(4x,A,i9)'    ,'Total rhs function calls   =',nfevals
    print '(4x,A,i9)'    ,'Num lin solver setup calls =',nlinsetups
    print '(4x,A,i9)'    ,'Num error test failures    =',netfails
    print '(4x,A,i9)'    ,'Last method order          =',qlast
    print '(4x,A,i9)'    ,'Next method order          =',qcur
    print '(4x,A,es12.5)','First internal step size   =',hinused
    print '(4x,A,es12.5)','Last internal step size    =',hlast
    print '(4x,A,es12.5)','Next internal step size    =',hcur
    print '(4x,A,es12.5)','Current internal time      =',tcur
    print '(4x,A,i9)'    ,'Num nonlinear solver iters =',nniters
    print '(4x,A,i9)'    ,'Num nonlinear solver fails =',nncfails
    print '(4x,A,i9)'    ,'Num Jacobian evaluations   =',njevals
    print *, ' '

    return

end subroutine CVodeStats
