program main
    use mod_prec
    use mod_implicit_solver, only: init, solve, destroy
    use mod_ode_analytic_boeuf
    use mpi
    implicit none

    real(pr) :: tstart, dtmin
    real(pr) :: rtol, atol
    integer :: ierr, i, j
    real(pr), dimension(:), allocatable :: y
    integer :: nloops = 100
    integer :: nstat  = 5
    real(pr) :: times, times_avg, times_var, times_stddev, stat_total
    real(pr), allocatable :: loop_times(:)

    call MPI_Init(ierr)

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
    call init(neq, rtol, atol, 0)

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
    call init( neq*nloops, rtol, atol, 0 )

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

    call MPI_Finalize(ierr)


end program main
