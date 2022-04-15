module mod_ode_analytic
use, intrinsic :: iso_c_binding

implicit none

integer(c_long), parameter :: neq = 8

! Parameters for each system type
real(8), dimension(:) :: Efield(2)                         !< Electric field at the current node
real(8), dimension(:) :: qom(10)                           !< Charge over mass for each species

! ODE_BENILOV
! CGB: possible to work with pointers here to avoid copies
real(8), dimension(:) :: species_mass_impl(2)              !< Mass of the species
real(8)               :: nu_ionization_impl                !< Ionization frequency
real(8)               :: int_electron_2_impl               !< Electron elastic collision frequency
real(8)               :: nu_ioniz_energy_electron_impl     !< Electron ionization energy frequency
real(8)               :: int_ion_impl                      !< Ion elastic collision frequency
real(8)               :: nu_ioniz_ion_impl                 !< Ion ionization frequency
real(8)               :: n_neutre_impl                     !< Neutral gas background density
real(8)               :: T_gas_impl                        !< Neutral gas background temperature

integer :: n_problems


contains

integer function RhsFn ( t, sunvec_y, sunvec_f, user_data ) result(ierr) bind(C, name="RhsFn_benilov")
    use fsundials_nvector_mod
    implicit none

    ! IN/OUT
    real(8),    value :: t
    type(N_Vector)     :: sunvec_y
    type(N_Vector)     :: sunvec_f
    type(c_ptr), value :: user_data

    ! LOCAL
    real(8), dimension(:), pointer :: yvec
    real(8), dimension(:), pointer :: fvec
    real(8), dimension(:,:), pointer :: y
    real(8), dimension(:,:), pointer :: f

    real(8) :: half, two, three, Boltzmann
    integer :: i

    ! Get data arrays from Sundials vector
    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)

    ! Pointer bounds remapping
    y(1:neq, 1:n_problems) => yvec
    f(1:neq, 1:n_problems) => fvec

    half = 0.5D0 ; two = 2.0D0 ; three = 3.0D0 ; Boltzmann = 1.380649D-23

    ! do i = 1, n_problems
    ! ! Fill RHS vector
    ! ! rho  ionization
    ! f(i) = n_neutre_impl * y(i) * nu_ionization_impl
    ! f(i+1) = n_neutre_impl * y(i) * nu_ioniz_ion_impl
    ! ! rhou    electric               + elastic collisions
    ! f(i+2) = qom(1) * y(i) * Efield(1) - n_neutre_impl * y(i+2) * int_electron_2_impl
    ! f(i+3) = qom(2) * y(i+1) * Efield(1) - n_neutre_impl * y(i+3) * int_ion_impl
    ! ! rhov    electric               + elastic collisions
    ! f(i+4) = qom(1) * y(i) * Efield(2) - n_neutre_impl * y(i+4) * int_electron_2_impl
    ! f(i+5) = qom(2) * y(i+1) * Efield(2) - n_neutre_impl * y(i+5) * int_ion_impl
    ! ! rhoE_e    electric - (excitation + ionization) - elastic collisions
    ! f(7, :) = qom(1) * (y(i+2) * Efield(1) + y(i+4) * Efield(2)) - n_neutre_impl * y(i) * nu_ioniz_energy_electron_impl &
    !         - two * species_mass_impl(1) / species_mass_impl(2) * n_neutre_impl * int_electron_2_impl &
    !             * ( y(i+6) - three * y(i) * Boltzmann * T_gas_impl / species_mass_impl(1) )
    ! ! rhoE_i   electric + ionization - elastic collisions
    ! f(i+7) = qom(2) * (y(i+3) * Efield(1) + y(i+5) * Efield(2)) &
    !         + three * half * Boltzmann * T_gas_impl * n_neutre_impl * y(i) * nu_ioniz_ion_impl &
    !         - n_neutre_impl * int_ion_impl * ( y(i+7) - three * half * y(i+1) * Boltzmann * T_gas_impl / species_mass_impl(2) )
    ! end do

    ! Fill RHS vector
    ! rho  ionization
    f(1, :) = n_neutre_impl * y(1, :) * nu_ionization_impl
    f(2, :) = n_neutre_impl * y(1, :) * nu_ioniz_ion_impl
    ! rhou    electric               + elastic collisions
    f(3, :) = qom(1) * y(1, :) * Efield(1) - n_neutre_impl * y(3, :) * int_electron_2_impl
    f(4, :) = qom(2) * y(2, :) * Efield(1) - n_neutre_impl * y(4, :) * int_ion_impl
    ! rhov    electric               + elastic collisions
    f(5, :) = qom(1) * y(1, :) * Efield(2) - n_neutre_impl * y(5, :) * int_electron_2_impl
    f(6, :) = qom(2) * y(2, :) * Efield(2) - n_neutre_impl * y(6, :) * int_ion_impl
    ! rhoE_e    electric - (excitation + ionization) - elastic collisions
    f(7, :) = qom(1) * (y(3, :) * Efield(1) + y(5, :) * Efield(2)) - n_neutre_impl * y(1, :) * nu_ioniz_energy_electron_impl &
            - two * species_mass_impl(1) / species_mass_impl(2) * n_neutre_impl * int_electron_2_impl &
                * ( y(7, :) - three * y(1, :) * Boltzmann * T_gas_impl / species_mass_impl(1) )
    ! rhoE_i   electric + ionization - elastic collisions
    f(8, :) = qom(2) * (y(4, :) * Efield(1) + y(6, :) * Efield(2)) &
            + three * half * Boltzmann * T_gas_impl * n_neutre_impl * y(1, :) * nu_ioniz_ion_impl &
            - n_neutre_impl * int_ion_impl * ( y(8, :) - three * half * y(2, :) * Boltzmann * T_gas_impl / species_mass_impl(2) )

    ierr = 0
    return

end function RhsFn

end module mod_ode_analytic


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
        use mod_ode_analytic, only: RhsFn

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

    subroutine solve ( y, tout )
        real(c_double), dimension(neq), intent(inout) :: y
        real(c_double), intent(in) :: tout
        real(c_double):: tcur(1)

        integer :: ierr
        real(8), pointer :: y_loc(:)

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
    use mod_ode_analytic
    use mpi
    implicit none

    real(8) :: tstart, dtmin
    real(8) :: rtol, atol
    integer :: ierr, i, j
    real(8), dimension(:), allocatable :: y
    integer :: nloops = 100
    integer :: nstat  = 5
    real(8) :: times, times_avg, times_var, times_stddev, stat_total
    real(8), allocatable :: loop_times(:)

    ! -----------------------------------------------------------------
    !    Serial
    ! -----------------------------------------------------------------

    allocate ( y(neq) )
    allocate ( loop_times(nloops) )

    tstart = 0.0D0
    dtmin = 1.0D-10

    rtol = 1.0D-10
    atol = 1.0D-50

    ! Set values
    Efield(1) = 1.0D5 ; Efield(2) = - 1.0D5
    qom(1) = - 1.0d-19 / 1.0D-32 ; qom(2) = - 1.0d-19 / 1.0D-27
    species_mass_impl(1) = 1.0D-32 ; species_mass_impl(2) = 1.0D-27
    nu_ionization_impl = 1.0D-14
    int_electron_2_impl = 1.0D-14
    nu_ioniz_energy_electron_impl = 1.0D-1
    int_ion_impl = 1.0D-15
    nu_ioniz_ion_impl = 1.0D-10
    n_neutre_impl = 9.64d20
    T_gas_impl = 300.0d0

    n_problems = 1

    print *, "init"
    call init(neq, rtol, atol)

    print *, "loop"
    stat_total = 0.0D0
    do j = 1, nstat
        do i = 1, nloops
            y(1) = 1.0D-15 ; y(2) = 1.0D-12
            y(3) = 1.5D-7  ; y(4) = 3.0D-16
            y(5) = 1.0D-20 ; y(6) = 1.0D-23
            y(7) = 1.0D-15 ; y(8) = 1.0D-15
            loop_times(i) = MPI_Wtime()
            call solve ( y, dtmin )
            loop_times(i) = MPI_Wtime() - loop_times(i)
        end do
        stat_total = stat_total + SUM(loop_times)
    end do

    times_avg = SUM(loop_times) / nloops
    times_var = 0.0D0
    do i = 1, nloops
        times_var = times_var + (loop_times(i) - times_avg)**2
    end do
    times_var = times_var / (nloops - 1)
    times_stddev = SQRT(times_var)

    write (*, '(A,ES15.8,A)') "time_average=", times_avg, " s"
    write (*, '(A,ES15.8,A)') "time_stddev=",  times_stddev, " s"
    write (*, '(A,ES15.8,A)') "time_best=", MINVAL(loop_times), " s"
    write (*, '(A,ES15.8,A)') "time_total=", stat_total / nstat, " s"


    call destroy()
    deallocate ( y )
    deallocate ( loop_times )

    ! ----------------------------------------------------------------------------
    !    Vectorized
    ! ----------------------------------------------------------------------------

    allocate ( y(neq*nloops) )
    allocate ( loop_times(1) )
    call init( neq*nloops, rtol, atol )

    n_problems = nloops
    stat_total = 0.0D0

    do j = 1, nstat
        do i = 1, nloops
            y(i) = 1.0D-15 ; y(i + 1) = 1.0D-12
            y(i + 2) = 1.5D-7  ; y(i + 3) = 3.0D-16
            y(i + 4) = 1.0D-20 ; y(i + 5) = 1.0D-23
            y(i + 6) = 1.0D-15 ; y(i + 7) = 1.0D-15
        end do
        print *, "vector solve"
        loop_times(1) = MPI_Wtime()
        call solve ( y, dtmin )
        loop_times(1) = MPI_Wtime() - loop_times(1)
        stat_total = stat_total + loop_times(1)
    end do

    write (*, '(A,ES15.8,A)') "time_total=", stat_total / nstat, " s"

    deallocate ( y )
    deallocate ( loop_times )
    call destroy ()


end program main
