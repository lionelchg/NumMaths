module mod_ode_analytic_boeuf
use, intrinsic :: iso_c_binding
use mod_prec

implicit none

integer(c_long), parameter :: neq = 8

! Parameters for each system type
real(pr), dimension(:) :: Efield(2)                         !< Electric field at the current node
real(pr), dimension(:) :: qom(10)                           !< Charge over mass for each species

! ODE_BENILOV
! CGB: possible to work with pointers here to avoid copies
real(pr), dimension(:) :: species_mass_impl(2)              !< Mass of the species
real(pr)               :: nu_ionization_impl                !< Ionization frequency
real(pr)               :: int_electron_2_impl               !< Electron elastic collision frequency
real(pr)               :: nu_ioniz_energy_electron_impl     !< Electron ionization energy frequency
real(pr)               :: int_ion_impl                      !< Ion elastic collision frequency
real(pr)               :: nu_ioniz_ion_impl                 !< Ion ionization frequency
real(pr)               :: n_neutre_impl                     !< Neutral gas background density
real(pr)               :: T_gas_impl                        !< Neutral gas background temperature

integer :: n_problems


contains

integer function RhsFn ( t, sunvec_y, sunvec_f, user_data ) result(ierr) bind(C, name="RhsFn_benilov")
    use fsundials_nvector_mod
    use mod_prec
    implicit none

    ! IN/OUT
    real(pr),    value :: t
    type(N_Vector)     :: sunvec_y
    type(N_Vector)     :: sunvec_f
    type(c_ptr), value :: user_data

    ! LOCAL
    real(pr), dimension(:), pointer :: yvec
    real(pr), dimension(:), pointer :: fvec
    real(pr), dimension(:,:), pointer :: y
    real(pr), dimension(:,:), pointer :: f

    real(pr) :: half, two, three, Boltzmann
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

! integer function JacFn_benilov(t, sunvec_y, sunvec_f, sunmat_J, &
!     user_data, tmp1, tmp2, tmp3) result(ierr) bind(C, name="JacFn_benilov")

!     use fsundials_nvector_mod
!     use fsunmatrix_dense_mod
!     use fsundials_matrix_mod
!     implicit none

!     ! IN/OUT
!     real(pr),    value :: t                     !< Current time
!     type(N_Vector)     :: sunvec_y              !< Sundials solution vector
!     type(N_Vector)     :: sunvec_f              !< Sundials RHS vector (~ ydot)
!     type(SUNMatrix)    :: sunmat_J              !< Jacobian SUNMatrix
!     type(c_ptr), value :: user_data             !< Pointer to user-defined data structure (not used)
!     type(N_Vector)     :: tmp1, tmp2, tmp3      !< Workspace N_Vectors

!     ! LOCAL
!     real(pr), dimension(:), pointer :: yvec     !< Pointer to data in N_Vector
!     real(pr), dimension(:), pointer :: fvec     !< Pointer to data in N_Vector
!     real(pr), dimension(:), pointer :: Jmat     !< Pointer to data in SUNMatrix

!     ! Get data arrays from Sundials vector and matrix
!     yvec => FN_VGetArrayPointer(sunvec_y)
!     Jmat => FSUNDenseMatrix_Data(sunmat_J)

!     ! Fill Jacobian matrix in column major order (each line is a column of the Jacobian)
!     ! Aliases: nu_ie -> nu_ioniz_electron; nu_ii -> nu_ioniz_ion; nu_ee -> nu_elastic_electron; nu_ei -> nu_elastic_ion
!     ! q -> qom; E -> Efield; Ie -> ioniz_rhoEe; Ii -> ioniz_rhoEi; Een -> elast_rhoEe_n; Ein -> elast_rhoEi_n
!     ! |     nu_ie     |    0    |    0   |    0   |    0   |    0   |   0  |   0  |
!     ! |     nu_ii     |    0    |    0   |    0   |    0   |    0   |   0  |   0  |
!     ! |     q1*E1     |    0    | -nu_ee |    0   |    0   |    0   |   0  |   0  |
!     ! |       0       |  q2*E1  |    0   | -nu_ei |    0   |    0   |   0  |   0  |
!     ! |     q1*E2     |    0    |    0   |    0   | -nu_ee |    0   |   0  |   0  |
!     ! |       0       |  q2*E2  |    0   |    0   |    0   | -nu_ei |   0  |   0  |
!     ! | -Ie + Ee1*Ee2 |    0    |  q1*E1 |    0   |  q1*E2 |    0   | -Ee1 |   0  |
!     ! |       Ii      | Ei1*Ei2 |    0   |  q2*E1 |    0   |  q2*E2 |   0  | -Ei1 |
!     Jmat = [ &
!         nu_ioniz_electron_impl, nu_ioniz_ion_impl, qom(1) * Efield(1), zero, qom(1) * Efield(2), zero, - ioniz_rhoEe + elast_rhoEe_1 * elast_rhoEe_2, ioniz_rhoEi, &
!         zero, zero, zero, qom(2) * Efield(1), zero, qom(2) * Efield(2), zero, elast_rhoEi_1 * elast_rhoEi_2,    &
!         zero, zero, - nu_elastic_electron_impl, zero, zero, zero, qom(1) * Efield(1), zero,                     &
!         zero, zero, zero, - nu_elastic_ion_impl, zero, zero, zero, qom(2) * Efield(1),                          &
!         zero, zero, zero, zero, - nu_elastic_electron_impl, zero, qom(1) * Efield(2), zero,                     &
!         zero, zero, zero, zero, zero, - nu_elastic_ion_impl, zero, qom(2) * Efield(2),                          &
!         zero, zero, zero, zero, zero, zero, - elast_rhoEe_1, zero,                                              &
!         zero, zero, zero, zero, zero, zero, zero, - elast_rhoEi_1                                               &
!     ]

!     ! call FSUNDenseMatrix_Print ( sunmat_J )

!     ! print *, "Jacobian"
!     ! write (*, '(8(ES15.8,2X))') nu_ioniz_electron_impl, nu_ioniz_ion_impl, qom(1) * Efield(1), zero, qom(1) * Efield(2), zero, - ioniz_rhoEe + elast_rhoEe_1 * elast_rhoEe_2, ioniz_rhoEi
!     ! write (*, '(8(ES15.8,2X))') zero, zero, zero, qom(2) * Efield(1), zero, qom(2) * Efield(2), zero, elast_rhoEi_1 * elast_rhoEi_2
!     ! write (*, '(8(ES15.8,2X))') zero, zero, - nu_elastic_electron_impl, zero, zero, zero, qom(1) * Efield(1), zero
!     ! write (*, '(8(ES15.8,2X))') zero, zero, zero, - nu_elastic_ion_impl, zero, zero, zero, qom(2) * Efield(1)
!     ! write (*, '(8(ES15.8,2X))') zero, zero, zero, zero, - nu_elastic_electron_impl, zero, qom(1) * Efield(2), zero
!     ! write (*, '(8(ES15.8,2X))') zero, zero, zero, zero, zero, - nu_elastic_ion_impl, zero, qom(2) * Efield(2)
!     ! write (*, '(8(ES15.8,2X))') zero, zero, zero, zero, zero, zero, - elast_rhoEe_1, zero
!     ! write (*, '(8(ES15.8,2X))') zero, zero, zero, zero, zero, zero, zero, - elast_rhoEi_1
!     ! print *, ""
!     ! print *, ""
!     ! print *, ""

!     ! Return success
!     ierr = 0
!     return

! end function JacFn


end module mod_ode_analytic_boeuf