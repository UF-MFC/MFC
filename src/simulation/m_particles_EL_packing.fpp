!>
!! @file m_particles_EL_packing.fpp
!! @brief Generates Euler-Lagrange particle initial conditions.

#:include 'macros.fpp'

!> @brief Deterministic rejection-sampling packer for Euler-Lagrange particles.
module m_particles_EL_packing

    use m_global_parameters
    use m_mpi_common
    use m_derived_types, only: ib_patch_parameters, particle_cloud_parameters
    use m_particle_cloud, only: s_generate_particle_clouds

    implicit none

    private

    public :: s_estimate_EL_particle_count
    public :: s_generate_EL_particle_input

contains

    !> Estimate number of particles from the requested packing volume fraction.
    impure subroutine s_estimate_EL_particle_count()

        real(wp) :: region_measure, particle_measure

        if (lag_params%packing_flag < 1) return

        call s_validate_EL_packing_params()

        particle_measure = f_mean_particle_measure()

        select case (lag_params%packing_flag)
        case (1)
            if (p == 0) then
                region_measure = lag_params%packing_length(1)*lag_params%packing_length(2)
            else
                region_measure = lag_params%packing_length(1)*lag_params%packing_length(2)*lag_params%packing_length(3)
            end if
        case (2)
            if (p == 0) then
                region_measure = pi*(lag_params%packing_shell_outer_radius**2._wp - lag_params%packing_shell_inner_radius**2._wp)
                region_measure = 0.5_wp*region_measure
            else
                region_measure = 4._wp/3._wp*pi*(lag_params%packing_shell_outer_radius**3._wp &
                                                 & - lag_params%packing_shell_inner_radius**3._wp)
                region_measure = 0.5_wp*region_measure
            end if
        case default
            call s_mpi_abort("Unsupported lag_params%packing_flag. Use 0=file, 1=box, or 2=hemisphere shell.")
        end select

        lag_params%nParticles_glb = max(1, int(lag_params%packing_volume_fraction*region_measure/particle_measure))

        if (proc_rank == 0) then
            print '(a,i0)', 'Estimated EL packed particles: ', lag_params%nParticles_glb
        end if

    end subroutine s_estimate_EL_particle_count

    !> Generate packed particle input rows: x y z vx vy vz radius dummy.
    impure subroutine s_generate_EL_particle_input(input_particles, n_generated)

        real(wp), allocatable, intent(out), dimension(:,:) :: input_particles
        integer, intent(out)                               :: n_generated
        integer                                            :: i, num_particle_clouds_save
        type(particle_cloud_parameters)                    :: particle_cloud_save
        type(ib_patch_parameters), allocatable             :: particle_cloud_ibs(:)
        real(wp)                                           :: radius, length_x, length_y, length_z

        allocate (input_particles(8, lag_params%nParticles_glb))

        radius = 0.5_wp*lag_params%packing_diameter_mean
        length_x = lag_params%packing_length(1)
        length_y = lag_params%packing_length(2)
        length_z = lag_params%packing_length(3)
        if (lag_params%packing_flag == 2) then
            length_x = 2._wp*lag_params%packing_shell_outer_radius
            length_y = 2._wp*lag_params%packing_shell_outer_radius
            length_z = 2._wp*lag_params%packing_shell_outer_radius
        end if

        num_particle_clouds_save = num_particle_clouds
        particle_cloud_save = particle_cloud(1)

        num_particle_clouds = 1
        particle_cloud(1)%x_centroid = lag_params%packing_centroid(1)
        particle_cloud(1)%y_centroid = lag_params%packing_centroid(2)
        particle_cloud(1)%z_centroid = lag_params%packing_centroid(3)
        particle_cloud(1)%length_x = length_x
        particle_cloud(1)%length_y = length_y
        particle_cloud(1)%length_z = length_z
        particle_cloud(1)%num_particles = lag_params%nParticles_glb
        particle_cloud(1)%radius = radius
        particle_cloud(1)%mass = 0._wp
        particle_cloud(1)%min_spacing = lag_params%packing_min_spacing
        particle_cloud(1)%shell_inner_radius = lag_params%packing_shell_inner_radius
        particle_cloud(1)%shell_outer_radius = lag_params%packing_shell_outer_radius
        particle_cloud(1)%moving_ibm = 0
        particle_cloud(1)%seed = lag_params%packing_seed
        particle_cloud(1)%packing_method = lag_params%packing_flag

        ! Map EL packing flags to particle_cloud packing methods.
        if (lag_params%packing_flag == 2) particle_cloud(1)%packing_method = 3

        call s_generate_particle_clouds(particle_cloud_ibs)
        n_generated = size(particle_cloud_ibs)

        do i = 1, n_generated
            input_particles(1, i) = particle_cloud_ibs(i)%x_centroid
            input_particles(2, i) = particle_cloud_ibs(i)%y_centroid
            input_particles(3, i) = particle_cloud_ibs(i)%z_centroid
            input_particles(4, i) = lag_params%packing_velocity(1)
            input_particles(5, i) = lag_params%packing_velocity(2)
            input_particles(6, i) = lag_params%packing_velocity(3)
            input_particles(7, i) = radius
            input_particles(8, i) = 0._wp
        end do

        deallocate (particle_cloud_ibs)

        num_particle_clouds = num_particle_clouds_save
        particle_cloud(1) = particle_cloud_save

    end subroutine s_generate_EL_particle_input

    subroutine s_validate_EL_packing_params()

        if (lag_params%packing_volume_fraction <= 0._wp) then
            call s_mpi_abort("lag_params%packing_volume_fraction must be positive for EL particle packing.")
        end if

        if (lag_params%packing_max_attempts <= 0) then
            call s_mpi_abort("lag_params%packing_max_attempts must be positive for EL particle packing.")
        end if

        if (lag_params%packing_size_distribution /= 0) then
            call s_mpi_abort("EL particle packing currently reuses Dan's constant-radius packing path; use packing_size_distribution=0.")
        end if

        if (lag_params%packing_diameter_mean <= 0._wp) then
            call s_mpi_abort("Constant EL particle packing requires positive lag_params%packing_diameter_mean.")
        end if

        if (lag_params%packing_flag == 2) then
            if (lag_params%packing_shell_inner_radius < 0._wp &
                & .or. lag_params%packing_shell_outer_radius <= lag_params%packing_shell_inner_radius) then
                call s_mpi_abort("Hemisphere shell packing requires 0 <= inner_radius < outer_radius.")
            end if
        end if

        if (lag_params%packing_diameter_mean + lag_params%packing_min_spacing <= 0._wp) then
            call s_mpi_abort("EL particle packing requires positive maximum diameter plus minimum spacing.")
        end if

    end subroutine s_validate_EL_packing_params

    function f_mean_particle_measure() result(measure)

        real(wp) :: measure, radius

        radius = 0.5_wp*lag_params%packing_diameter_mean
        if (p == 0) then
            measure = pi*radius**2._wp
        else
            measure = 4._wp/3._wp*pi*radius**3._wp
        end if

    end function f_mean_particle_measure

end module m_particles_EL_packing
