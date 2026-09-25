!>
!! @file
!! @brief Contains module m_reactive_burn

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Bounded pressure burn and Garno ignition-and-growth source for the multi-fluid model. Garno uses air (fluid 1), unreacted
!! explosive (fluid 2), and products (fluid 3). The last two share an EOS; their combined partial density and volume fraction
!! represent one material.
module m_reactive_burn

    use m_global_parameters
    use m_variables_conversion, only: s_compute_mixture_coefficients
    use m_eos, only: s_phase_temperature, s_phase_internal_energy, f_pressure
    use m_constants, only: eos_jwl
    use m_ibm, only: ib_markers
    use m_mpi_proxy, only: s_mpi_abort

    implicit none

    private; public :: s_reactive_burn_substep, s_program_burn_step

contains

    !> Pressure-law burn frequency, with an optional reactant-temperature factor.
    subroutine s_burn_frequency(pres, alpha_rho_react, alpha_react, frequency)

        $:GPU_ROUTINE(function_name='s_burn_frequency', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: pres, alpha_rho_react, alpha_react
        real(wp), intent(out) :: frequency
        real(wp)              :: drive, T_r

        ! Pressure-driven burn fires only behind the shock (p > rburn%pign).
        drive = (pres - rburn%pign)/rburn%pref
        if (drive > 0._wp .and. alpha_react > 0._wp .and. alpha_rho_react > 0._wp) then
            frequency = rburn%k*drive**rburn%n
            ! Optional Arrhenius dependence on the reactant's phasic temperature, from the reactant's own
            ! EOS. rburn%ta = 0, the default, leaves the pure pressure-driven rate unchanged.
            if (rburn%ta > 0._wp) then
                call s_phase_temperature(alpha_rho_react/max(alpha_react, sgm_eps), pres, 1, T_r)
                frequency = frequency*exp(-rburn%ta/T_r)
            end if
        else
            frequency = 0._wp
        end if

    end subroutine s_burn_frequency

    !> Garno's RI + RG. Since dY_R/dt = -(RI + RG), this is a mass-fraction rate, not a frequency.
    subroutine s_garno_rate(rho, alpha_rho_react, rate)

        $:GPU_ROUTINE(function_name='s_garno_rate', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: rho, alpha_rho_react
        real(wp), intent(out) :: rate
        real(wp)              :: y_react, density_ratio

        rate = 0._wp
        if (rho <= sgm_eps .or. alpha_rho_react <= 0._wp) return
        y_react = min(max(alpha_rho_react/rho, 0._wp), 1._wp)
        density_ratio = rho/rburn%rho0
        rate = rburn%ki*y_react**rburn%m1*abs(1._wp - density_ratio)**rburn%m2 + rburn%kg*y_react**rburn%n1*(1._wp - y_react) &
                                              & **rburn%n2*density_ratio**rburn%n3

    end subroutine s_garno_rate

    !> Bounded update; the common Garno exponents reduce to an exactly solvable logistic equation.
    subroutine s_garno_survival(rho, alpha_rho_react, dtime, survival)

        $:GPU_ROUTINE(function_name='s_garno_survival', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: rho, alpha_rho_react, dtime
        real(wp), intent(out) :: survival
        real(wp)              :: y_react, y_mid, rate, a, b, c, decay

        survival = 1._wp
        if (rho <= sgm_eps .or. alpha_rho_react <= 0._wp) return
        y_react = min(max(alpha_rho_react/rho, 0._wp), 1._wp)

        if (rburn%m1 == 1._wp .and. rburn%n1 == 1._wp .and. rburn%n2 == 1._wp) then
            a = rburn%ki*abs(1._wp - rho/rburn%rho0)**rburn%m2
            b = rburn%kg*(rho/rburn%rho0)**rburn%n3
            c = a + b
            if (c <= 0._wp .or. (a == 0._wp .and. y_react >= 1._wp)) return
            decay = exp(-c*dtime)
            survival = min(decay/max(1._wp - (b/c)*y_react*(1._wp - decay), sgm_eps), 1._wp)
        else
            call s_garno_rate(rho, alpha_rho_react, rate)
            if (rate <= 0._wp) return
            y_mid = y_react*exp(-0.5_wp*rate*dtime/max(y_react, sgm_eps))
            if (y_mid <= sgm_eps) then
                survival = 0._wp
            else
                call s_garno_rate(rho, rho*y_mid, rate)
                survival = exp(-rate*dtime/y_mid)
            end if
        end if

    end subroutine s_garno_survival

    !> Integrate the local burn with a bounded source update.
    !! @param q_cons_vf  Conserved variables, updated in place
    !! @param dtime      Time step to integrate across
    !! @param bounds     Interior cell bounds
    subroutine s_reactive_burn_substep(q_cons_vf, dtime, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: dtime
        type(int_bounds_info), dimension(1:3), intent(in)      :: bounds
        integer                                                :: x, y, z, i, sub, nsub
        real(wp)                                               :: rho, pres, frequency, survival, denergy
        real(wp)                                               :: dt_sub, e_int, gamma_mix, pi_inf_mix, qv_mix
        real(wp)                                               :: rho_mix, dlambda, dmass, e_phase
        real(wp)                                               :: alpha_rho_react, alpha_react, alpha_phase, alpha_rho_phase

        ! num_fluids_max also covers the literal product index 3 under case optimization.
        real(wp), dimension(num_fluids_max) :: alpha_rho, alpha

        nsub = max(1, rburn%substeps)
        dt_sub = dtime/real(nsub, wp)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho, alpha, rho, pres, frequency, survival, denergy, e_int, gamma_mix, &
                            & pi_inf_mix, qv_mix, rho_mix, dlambda, dmass, e_phase, alpha_rho_react, alpha_react, alpha_phase, &
                            & alpha_rho_phase, i, sub]', copyin='[bounds, dt_sub, nsub]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(x, y, z)
                        alpha(i) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(x, y, z)
                    end do
                    rho = alpha_rho(1) + alpha_rho(2)
                    if (rburn%model == 1) rho = rho + alpha_rho(3)
                    if (rho <= sgm_eps) cycle

                    if (rburn%model == 0 .or. model_eqns == model_eqns_6eq) then
                        e_int = q_cons_vf(eqn_idx%E)%sf(x, y, z)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            e_int = e_int - 0.5_wp*q_cons_vf(i)%sf(x, y, z)**2/rho
                        end do
                    end if

                    $:GPU_LOOP(parallelism='[seq]')
                    do sub = 1, nsub
                        if (rburn%model == 0) then
                            call s_compute_mixture_coefficients(alpha_rho, alpha, rho_mix, gamma_mix, pi_inf_mix, qv_mix)
                            pres = f_pressure(e_int, gamma_mix, pi_inf_mix, qv_mix)
                            alpha_rho_react = alpha_rho(1)
                            alpha_react = alpha(1)
                            call s_burn_frequency(pres, alpha_rho_react, alpha_react, frequency)
                            if (frequency <= 0._wp) exit
                            survival = exp(-frequency*dt_sub)
                            dlambda = alpha(1)*(1._wp - survival)
                            dmass = alpha_rho(1)*(1._wp - survival)
                            alpha(1) = alpha(1)*survival
                            alpha(2) = alpha(2) + dlambda
                            alpha_rho(1) = alpha_rho(1)*survival
                            alpha_rho(2) = alpha_rho(2) + dmass
                        else
                            call s_garno_survival(rho, alpha_rho(2), dt_sub, survival)
                            if (survival >= 1._wp) exit
                            dlambda = alpha(2)*(1._wp - survival)
                            dmass = alpha_rho(2)*(1._wp - survival)
                            alpha(2) = alpha(2)*survival
                            alpha(3) = alpha(3) + dlambda
                            alpha_rho(2) = alpha_rho(2)*survival
                            alpha_rho(3) = alpha_rho(3) + dmass
                            denergy = rburn%q*dmass
                            q_cons_vf(eqn_idx%E)%sf(x, y, z) = q_cons_vf(eqn_idx%E)%sf(x, y, z) + denergy
                            if (model_eqns == model_eqns_6eq) e_int = e_int + denergy
                        end if
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(x, y, z) = alpha_rho(i)
                        q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(x, y, z) = alpha(i)
                    end do
                    if (model_eqns == model_eqns_6eq) then
                        call s_compute_mixture_coefficients(alpha_rho, alpha, rho_mix, gamma_mix, pi_inf_mix, qv_mix)
                        pres = f_pressure(e_int, gamma_mix, pi_inf_mix, qv_mix)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_phase = alpha(i)
                            alpha_rho_phase = alpha_rho(i)
                            call s_phase_internal_energy(pres, alpha_phase, alpha_rho_phase, i, e_phase)
                            q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(x, y, z) = e_phase
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_reactive_burn_substep

    !> Deposit the energy swept by the prescribed front during this flow step.
    subroutine s_program_burn_step(q_cons_vf, t_begin, dtime)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: t_begin, dtime
        real(wp)                                               :: r_start, r_end, r_cell, f_start, f_end, qdet, denergy
        real(wp)                                               :: x_det, y_det, z_det, width
        integer                                                :: i, jwl_phase, x, y, z, marker

        if (.not. prog_burn) return

        jwl_phase = 0
        do i = 1, num_fluids
            if (fluid_pp(i)%eos == eos_jwl) then
                if (jwl_phase /= 0) call s_mpi_abort('prog_burn requires exactly one JWL fluid')
                jwl_phase = i
            end if
        end do
        if (jwl_phase == 0) call s_mpi_abort('prog_burn requires one JWL fluid')
        qdet = fluid_pp(jwl_phase)%jwl_Q
        if (qdet <= 0._wp .or. pb_D_cj <= 0._wp .or. pb_width <= 0._wp) then
            call s_mpi_abort('prog_burn requires positive jwl_Q, pb_D_cj, and pb_width')
        end if

        r_start = pb_D_cj*max(t_begin - pb_t_det, 0._wp)
        r_end = pb_D_cj*max(t_begin + dtime - pb_t_det, 0._wp)
        if (r_end <= r_start) return
        x_det = pb_x_det
        y_det = pb_y_det
        z_det = pb_z_det
        width = pb_width

        $:GPU_PARALLEL_LOOP(collapse=3, private='[x, y, z, marker, r_cell, f_start, f_end, denergy]', copyin='[r_start, r_end, &
                            & qdet, jwl_phase, x_det, y_det, z_det, width]')
        do z = 0, p
            do y = 0, n
                do x = 0, m
                    marker = 0
                    if (ib) marker = ib_markers%sf(x, y, z)
                    if (marker == 0) then
                        r_cell = (x_cc(x) - x_det)**2
                        if (n > 0) r_cell = r_cell + (y_cc(y) - y_det)**2
                        if (p > 0) r_cell = r_cell + (z_cc(z) - z_det)**2
                        r_cell = sqrt(r_cell)
                        f_start = min(max((r_start - r_cell)/width, 0._wp), 1._wp)
                        f_end = min(max((r_end - r_cell)/width, 0._wp), 1._wp)
                        denergy = q_cons_vf(eqn_idx%cont%beg + jwl_phase - 1)%sf(x, y, z)*qdet*(f_end - f_start)
                        q_cons_vf(eqn_idx%E)%sf(x, y, z) = q_cons_vf(eqn_idx%E)%sf(x, y, z) + denergy
                        if (model_eqns == model_eqns_6eq) then
                            q_cons_vf(eqn_idx%int_en%beg + jwl_phase - 1)%sf(x, y, &
                                      & z) = q_cons_vf(eqn_idx%int_en%beg + jwl_phase - 1)%sf(x, y, z) + denergy
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_program_burn_step

end module m_reactive_burn
