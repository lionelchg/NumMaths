program main
    use mod_prec
    use mod_implicit_solver, only: init, solve, destroy
    use mod_ode_analytic_brusselator
    use mpi
    implicit none

    real(pr) :: tstart, dtmin
    real(pr) :: rtol, atol
    integer :: ierr, i, j, loc
    real(pr), dimension(:), allocatable :: y
    integer :: nloops
    integer :: nstat  = 5
    real(pr) :: times, times_avg, times_var, times_stddev, stat_total
    real(pr), allocatable :: loop_times(:)
    character(len=100) :: cli
    logical :: use_jac

    call MPI_Init(ierr)

    write (*,'(A)') " ----------------------------------------------------------------------------"
    write (*,'(A)') "                            Input parameters"
    write (*,'(A)') " ----------------------------------------------------------------------------"

    call get_command_argument (1, cli)
    read (cli, *) nloops
    write (*, '(A,I0)') "nloops=", nloops
    call get_command_argument (2, cli)
    if (cli == "true") then
        use_jac = .TRUE.
    else
        use_jac = .FALSE.
    end if
    write (*, '(A,L)') "use_jac=", use_jac

    write (*,'(A)') " ----------------------------------------------------------------------------"
    write (*,'(A)') "                                 Serial"
    write (*,'(A)') " ----------------------------------------------------------------------------"

    allocate ( y(neq) )
    allocate ( loop_times(nloops) )

    tstart = 0.0D0
    dtmin = 1.0D0

    rtol = 1.0D-10
    atol = 1.0D-20

    ! Set values
    n_problems = 1

    if ( use_jac ) then
        call init ( neq, rtol, atol, 2 )
    else
        call init ( neq, rtol, atol, 1 )
    end if

    ! write (*,'(A)') "loop"
    stat_total = 0.0D0
    do j = 1, nstat
        do i = 1, nloops
            y(1) = 3.9D0 ; y(2) = 1.1D0 ; y(3) = 2.8D0

            loop_times(i) = MPI_Wtime()
            call solve ( y, dtmin )
            loop_times(i) = MPI_Wtime() - loop_times(i)
        end do
        stat_total = stat_total + SUM(loop_times)
    end do

    times_avg = SUM(loop_times) / nloops
    times_var = 0.0D0
    if ( nloops > 1 ) then
        do i = 1, nloops
            times_var = times_var + (loop_times(i) - times_avg)**2
        end do
        times_var = times_var / (nloops - 1)
        times_stddev = SQRT(times_var)
    else
        times_stddev = 0.0D0
    end if

    write (*, '(A,ES15.8,A)') "time_average=", times_avg, " s"
    write (*, '(A,ES15.8,A)') "time_stddev=",  times_stddev, " s"
    write (*, '(A,ES15.8,A)') "time_best=", MINVAL(loop_times), " s"
    write (*, '(A,ES15.8,A)') "time_total=", stat_total, " s"

    call destroy()
    deallocate ( y )
    deallocate ( loop_times )

    write (*,'(A)') ""
    write (*,'(A)') " ----------------------------------------------------------------------------"
    write (*,'(A)') "                               Vectorized"
    write (*,'(A)') " ----------------------------------------------------------------------------"

    n_problems = nloops
    allocate ( y(neq*n_problems) )
    allocate ( loop_times(1) )
    if (use_jac) then
        call init ( neq*n_problems, rtol, atol, 2 )
    else
        call init ( neq*n_problems, rtol, atol, 1 )
    end if


    stat_total = 0.0D0
    ! write (*,'(A)') "vector solve"
    do j = 1, nstat
        do i = 1, n_problems
            loc = (i-1) * neq
            y(loc + 1) = 3.9D0 ; y(loc + 2) = 1.1D0 ; y(loc + 3) = 2.8D0
        end do

        loop_times(1) = MPI_Wtime()
        call solve ( y, dtmin )
        loop_times(1) = MPI_Wtime() - loop_times(1)

        stat_total = stat_total + loop_times(1)
    end do

    write (*, '(A,ES15.8,A)') "time_total=", stat_total, " s"
    write (*,'(A)') ""

    deallocate ( y )
    deallocate ( loop_times )
    call destroy ()

    call MPI_Finalize(ierr)


end program main
