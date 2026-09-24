!>
!! @file
!! @brief Contains module m_lso_pp_filter

#:include 'macros.fpp'

!> @brief Post_process-side additional Gaussian filter.
!!
!! Applies 9-point FIR passes with total width sigma_2 = sqrt(sigma_target^2 - sigma_in^2)
!! using the lso_pp_a_* weights from the Python BCD design. Wide targets use the stage-2
!! lso_pp2_a_* weights as a second equal-variance cascade.
!!
!! IB-aware normalization uses the simulation-filtered mask lso_mask_<t>.dat, so that
!! filter2(w*qhat1)/filter2(w) = filter2(filter1(m*q))/filter2(filter1(m)).
module m_lso_pp_filter

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_constants
    use m_variables_conversion, only: gammas, pi_infs, qvs
    use m_eos, only: f_isentrope_pressure

    implicit none

    private

    public :: s_initialize_lso_pp_filter_module, s_finalize_lso_pp_filter_module, s_apply_lso_pp_filter, &
        & s_compute_lso_pp_stat_fields, s_filter_lso_pp_stat_fields, q_lso_pp_stat_vf, s_apply_lso_pp_filter_masked, &
        & s_compute_lso_closure_fields, f_lso_n_closure, q_lso_pp_w_vf

    ! Floor on the normalized-convolution denominator filter(w).
    real(wp), parameter :: lso_w_floor = 1.0e-3_wp
    integer, parameter  :: lso_pp_filter_radius = 4

    ! Scratch buffer for one directional pass.
    real(wp), allocatable :: lso_pp_tmp(:,:,:)

    ! Simulation-filtered statistical products, optionally widened in post-process.
    type(scalar_field), allocatable :: q_lso_pp_stat_vf(:)

    ! Post-filtered normalization weight w2 = filter2(w1) (interior only), used by the closure pass; 1 without a mask.
    type(scalar_field) :: q_lso_pp_w_vf(1:1)

contains

    impure subroutine s_initialize_lso_pp_filter_module()

        integer :: i

        if (lso_pp_filter) allocate (lso_pp_tmp(0:m,0:n,0:p))

        if (lso_pp_filter .or. lso_closure_wrt) then
            allocate (q_lso_pp_w_vf(1)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end, &
                      & idwbuff(3)%beg:idwbuff(3)%end))
            q_lso_pp_w_vf(1)%sf = 1._stp
        end if

        if (lso_stat_wrt .and. n_lso_stat > 0) then
            allocate (q_lso_pp_stat_vf(1:n_lso_stat))
            do i = 1, n_lso_stat
                allocate (q_lso_pp_stat_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end, &
                          & idwbuff(3)%beg:idwbuff(3)%end))
            end do
        end if

    end subroutine s_initialize_lso_pp_filter_module

    impure subroutine s_finalize_lso_pp_filter_module()

        integer :: i

        if (allocated(lso_pp_tmp)) deallocate (lso_pp_tmp)
        if (associated(q_lso_pp_w_vf(1)%sf)) deallocate (q_lso_pp_w_vf(1)%sf)

        if (allocated(q_lso_pp_stat_vf)) then
            do i = 1, n_lso_stat
                if (associated(q_lso_pp_stat_vf(i)%sf)) deallocate (q_lso_pp_stat_vf(i)%sf)
            end do
            deallocate (q_lso_pp_stat_vf)
        end if

    end subroutine s_finalize_lso_pp_filter_module

    !> Apply the post_process LSO filter to an interior-only field in place.
    impure subroutine s_apply_lso_pp_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        type(scalar_field), allocatable   :: q_work_vf(:)
        integer                           :: i, j, k, l, work_buff_size

        work_buff_size = max(buff_size, lso_pp_filter_radius)
        allocate (q_work_vf(1:size(q_cons_vf)))
        do i = 1, size(q_cons_vf)
            allocate (q_work_vf(i)%sf(-work_buff_size:m + work_buff_size,-work_buff_size:n + work_buff_size, &
                      & -work_buff_size:p + work_buff_size))
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_work_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        call s_apply_lso_pp_filter_ghosted(q_work_vf)
        if (lso_pp2_n_passes_x > 0 .or. lso_pp2_n_passes_y > 0 .or. lso_pp2_n_passes_z > 0) then
            call s_apply_lso_pp_filter_ghosted(q_work_vf, 2)
        end if

        do i = 1, size(q_cons_vf)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = q_work_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
            deallocate (q_work_vf(i)%sf)
        end do
        deallocate (q_work_vf)
        call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        if (n > 0) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
        if (p > 0) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)

    end subroutine s_apply_lso_pp_filter

    !> Apply the post_process LSO filter to a ghosted work field.
    impure subroutine s_apply_lso_pp_filter_ghosted(q_cons_vf, stage)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer, intent(in), optional     :: stage
        integer                           :: i, ipass, j, k, l, nv
        integer                           :: stage_number, npasses_x, npasses_y, npasses_z
        real(wp)                          :: a_x(5, lso_max_passes), a_y(5, lso_max_passes), a_z(5, lso_max_passes)
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_cons_vf)
        stage_number = 1
        if (present(stage)) stage_number = stage
        if (stage_number == 2) then
            npasses_x = lso_pp2_n_passes_x
            npasses_y = lso_pp2_n_passes_y
            npasses_z = lso_pp2_n_passes_z
            a_x = lso_pp2_a_x
            a_y = lso_pp2_a_y
            a_z = lso_pp2_a_z
        else
            npasses_x = lso_pp_n_passes_x
            npasses_y = lso_pp_n_passes_y
            npasses_z = lso_pp_n_passes_z
            a_x = lso_pp_a_x
            a_y = lso_pp_a_y
            a_z = lso_pp_a_z
        end if

        ! x-direction

        call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        do ipass = 1, npasses_x
            c0 = a_x(1, ipass)
            c1 = a_x(2, ipass)
            c2 = a_x(3, ipass)
            c3 = a_x(4, ipass)
            c4 = a_x(5, ipass)
            do i = 1, nv
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j - 1, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_cons_vf(i)%sf(j - 2, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_cons_vf(i)%sf(j - 3, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_cons_vf(i)%sf(j - 4, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            if (ipass < npasses_x) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        end do

        ! y-direction (2D/3D)
        if (n > 0) then
            call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
            do ipass = 1, npasses_y
                c0 = a_y(1, ipass)
                c1 = a_y(2, ipass)
                c2 = a_y(3, ipass)
                c3 = a_y(4, ipass)
                c4 = a_y(5, ipass)
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k - 1, &
                                           & l), wp) + real(q_cons_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_cons_vf(i)%sf(j, &
                                           & k - 2, l), wp) + real(q_cons_vf(i)%sf(j, k + 2, l), &
                                           & wp)) + c3*(real(q_cons_vf(i)%sf(j, k - 3, l), wp) + real(q_cons_vf(i)%sf(j, k + 3, &
                                           & l), wp)) + c4*(real(q_cons_vf(i)%sf(j, k - 4, l), wp) + real(q_cons_vf(i)%sf(j, &
                                           & k + 4, l), wp))
                            end do
                        end do
                    end do
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < npasses_y) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
            end do
        end if

        ! z-direction (3D)
        if (p > 0) then
            call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)
            do ipass = 1, npasses_z
                c0 = a_z(1, ipass)
                c1 = a_z(2, ipass)
                c2 = a_z(3, ipass)
                c3 = a_z(4, ipass)
                c4 = a_z(5, ipass)
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k, &
                                           & l - 1), wp) + real(q_cons_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_cons_vf(i)%sf(j, &
                                           & k, l - 2), wp) + real(q_cons_vf(i)%sf(j, k, l + 2), &
                                           & wp)) + c3*(real(q_cons_vf(i)%sf(j, k, l - 3), wp) + real(q_cons_vf(i)%sf(j, k, &
                                           & l + 3), wp)) + c4*(real(q_cons_vf(i)%sf(j, k, l - 4), wp) + real(q_cons_vf(i)%sf(j, &
                                           & k, l + 4), wp))
                            end do
                        end do
                    end do
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < npasses_z) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)
            end do
        end if

    end subroutine s_apply_lso_pp_filter_ghosted

    !> Mask-normalized filtering: q <- filter(w*q)/filter(w), applied in place (the solid interior takes the fluid average).
    impure subroutine s_apply_lso_pp_filter_masked(q_cons_vf, w_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:), w_vf(1:1)
        integer                           :: i, j, k, l, nv

        nv = size(q_cons_vf)
        do i = 1, nv
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), wp)*real(w_vf(1)%sf(j, k, l), wp), stp)
                    end do
                end do
            end do
        end do

        call s_apply_lso_pp_filter(q_cons_vf)
        call s_apply_lso_pp_filter(w_vf)

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    q_lso_pp_w_vf(1)%sf(j, k, l) = w_vf(1)%sf(j, k, l)
                end do
            end do
        end do

        do i = 1, nv
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), wp)/max(real(w_vf(1)%sf(j, k, l), wp), &
                                  & lso_w_floor), stp)
                    end do
                end do
            end do
        end do

    end subroutine s_apply_lso_pp_filter_masked

    !> Number of closure output fields for the current dimensionality: R_sg (nt) | Q_T (nd) | E_ku (nd) | W_tau_u (nd) | R_mu_sg
    !! (nt) | R_lam_sg (nd) | T_tilde (1) | u_favre (nd), with nt = 1/3/6 symmetric tensor components.
    pure integer function f_lso_n_closure() result(n_cls)

        integer :: nt

        nt = num_dims*(num_dims + 1)/2
        n_cls = 2*nt + 5*num_dims + 1

    end function f_lso_n_closure

    !> Euler-Lagrange closure fields from the LSO stat products, with F[.] the phase-weighted filter, rho_b = F[rho], u~ = F[rho
    !! u]/rho_b and T~ from the central EOS energy reconstruction, F[E] = E_lso*w: R_sg,ij = F[rho u_i u_j] - F[rho u_i] F[rho
    !! u_j]/rho_b Q_T,i = gamma Cv (F[rho u_i T] - T~ F[rho u_i]) E_ku,i = (F[rho u_i |u|^2] - F[rho u_i] F[rho |u|^2]/rho_b)/2
    !! W_tau_u,i = (F[rho (tau.u)_i] - (tau_b . F[rho u])_i)/rho_b R_mu_sg = tau_b - mu [grad(u~) + grad(u~)^T - (2/3) div(u~) I]
    !! R_lam_sg = q_b + lambda grad(T~), with the intensive averages tau_b = F[g tau]/w and q_b = F[g q]/w. Gradients are centred
    !! differences on the output grid with ghost-refreshed Favre fields; layout per f_lso_n_closure().
    impure subroutine s_compute_lso_closure_fields(q_stat_vf, q_cons_vf, w_vf, q_cls_vf)

        type(scalar_field), intent(in)    :: q_stat_vf(:), q_cons_vf(:), w_vf(1:1)
        type(scalar_field), intent(inout) :: q_cls_vf(:)
        type(scalar_field)                :: uT_vf(1:4)
        integer                           :: nd, nt, i, c, j, k, l, a, b, work_buff_size
        integer                           :: ti(6), tj(6)
        integer                           :: i_rsg, i_qt, i_eku, i_wtu, i_rmu, i_rlam, i_tt, i_uf
        real(wp)                          :: Cv, gcv, rho_b, F_E, wloc, dsp(3), stiffness
        real(wp)                          :: rhou(3), tdotru(3), dudx(3, 3), dT(3), div_u, tau_res, tb(3, 3)

        nd = num_dims
        nt = nd*(nd + 1)/2
        ti(1:6) = (/1, 1, 1, 2, 2, 3/)
        tj(1:6) = (/1, 2, 3, 2, 3, 3/)
        if (nd == 2) then
            ti(1:3) = (/1, 1, 2/); tj(1:3) = (/1, 2, 2/)
        end if
        i_rsg = 0; i_qt = i_rsg + nt; i_eku = i_qt + nd; i_wtu = i_eku + nd
        i_rmu = i_wtu + nd; i_rlam = i_rmu + nt; i_tt = i_rlam + nd; i_uf = i_tt + 1

        Cv = lso_R_gas*gammas(1)
        gcv = (1._wp + 1._wp/gammas(1))*Cv
        stiffness = f_isentrope_pressure(pi_infs(1), gammas(1))

        ! Favre velocity and temperature with ghost extents for MPI exchange and gradient stencils.
        work_buff_size = max(buff_size, lso_pp_filter_radius)
        do i = 1, nd + 1
            allocate (uT_vf(i)%sf(-work_buff_size:m + work_buff_size,-work_buff_size:n + work_buff_size, &
                      & -work_buff_size:p + work_buff_size))
        end do
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho_b = max(real(q_stat_vf(lso_stat_rho_beg)%sf(j, k, l), wp), sgm_eps)
                    do a = 1, nd
                        uT_vf(a)%sf(j, k, l) = real(real(q_stat_vf(lso_stat_rhou_beg + a - 1)%sf(j, k, l), wp)/rho_b, stp)
                    end do
                    wloc = max(real(w_vf(1)%sf(j, k, l), wp), lso_w_floor)
                    F_E = real(q_cons_vf(eqn_idx%E)%sf(j, k, l), wp)*wloc
                    uT_vf(nd + 1)%sf(j, k, l) = real((F_E - 0.5_wp*real(q_stat_vf(lso_stat_rhoke_beg)%sf(j, k, l), &
                          & wp) - stiffness*wloc - qvs(1)*rho_b)/(Cv*rho_b), stp)
                end do
            end do
        end do
        call s_lso_pp_filter_ghost_refresh(uT_vf(1:nd + 1), 1)
        if (n > 0) call s_lso_pp_filter_ghost_refresh(uT_vf(1:nd + 1), 2)
        if (p > 0) call s_lso_pp_filter_ghost_refresh(uT_vf(1:nd + 1), 3)

        ! Quasi-uniform output-grid spacing (mean; exact for integer decimation).
        dsp = 1._wp
        dsp(1) = (x_cc(m) - x_cc(0))/real(max(m, 1), wp)
        if (n > 0) dsp(2) = (y_cc(n) - y_cc(0))/real(n, wp)
        if (p > 0) dsp(3) = (z_cc(p) - z_cc(0))/real(p, wp)

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho_b = max(real(q_stat_vf(lso_stat_rho_beg)%sf(j, k, l), wp), sgm_eps)
                    rhou = 0._wp
                    do a = 1, nd
                        rhou(a) = real(q_stat_vf(lso_stat_rhou_beg + a - 1)%sf(j, k, l), wp)
                    end do
                    ! The viscous products are phase-weighted (F[g tau], F[g q]); divide by the
                    ! filtered gas mask to recover the intensive gas-conditional averages.
                    wloc = max(real(w_vf(1)%sf(j, k, l), wp), lso_w_floor)
                    tb = 0._wp
                    do c = 1, nt
                        tb(ti(c), tj(c)) = real(q_stat_vf(lso_stat_tau_beg + c - 1)%sf(j, k, l), wp)/wloc
                        tb(tj(c), ti(c)) = tb(ti(c), tj(c))
                    end do

                    ! gradients of the Favre fields (centred; ghosts filled above)
                    dudx = 0._wp; dT = 0._wp
                    do a = 1, nd
                        dudx(a, 1) = (real(uT_vf(a)%sf(j + 1, k, l), wp) - real(uT_vf(a)%sf(j - 1, k, l), wp))/(2._wp*dsp(1))
                        if (n > 0) dudx(a, 2) = (real(uT_vf(a)%sf(j, k + 1, l), wp) - real(uT_vf(a)%sf(j, k - 1, l), &
                            & wp))/(2._wp*dsp(2))
                        if (p > 0) dudx(a, 3) = (real(uT_vf(a)%sf(j, k, l + 1), wp) - real(uT_vf(a)%sf(j, k, l - 1), &
                            & wp))/(2._wp*dsp(3))
                    end do
                    dT(1) = (real(uT_vf(nd + 1)%sf(j + 1, k, l), wp) - real(uT_vf(nd + 1)%sf(j - 1, k, l), wp))/(2._wp*dsp(1))
                    if (n > 0) dT(2) = (real(uT_vf(nd + 1)%sf(j, k + 1, l), wp) - real(uT_vf(nd + 1)%sf(j, k - 1, l), &
                        & wp))/(2._wp*dsp(2))
                    if (p > 0) dT(3) = (real(uT_vf(nd + 1)%sf(j, k, l + 1), wp) - real(uT_vf(nd + 1)%sf(j, k, l - 1), &
                        & wp))/(2._wp*dsp(3))
                    div_u = 0._wp
                    do a = 1, nd
                        div_u = div_u + dudx(a, a)
                    end do

                    ! R_sg
                    do c = 1, nt
                        q_cls_vf(i_rsg + c)%sf(j, k, l) = real(real(q_stat_vf(lso_stat_rhouu_beg + c - 1)%sf(j, k, l), &
                                 & wp) - rhou(ti(c))*rhou(tj(c))/rho_b, stp)
                    end do
                    ! Q_T, E_ku
                    do a = 1, nd
                        q_cls_vf(i_qt + a)%sf(j, k, l) = real(gcv*(real(q_stat_vf(lso_stat_rhouT_beg + a - 1)%sf(j, k, l), &
                                 & wp) - real(uT_vf(nd + 1)%sf(j, k, l), wp)*rhou(a)), stp)
                        q_cls_vf(i_eku + a)%sf(j, k, l) = real(0.5_wp*(real(q_stat_vf(lso_stat_rhouke_beg + a - 1)%sf(j, k, l), &
                                 & wp) - rhou(a)/rho_b*real(q_stat_vf(lso_stat_rhoke_beg)%sf(j, k, l), wp)), stp)
                    end do
                    ! W_tau_u
                    tdotru = 0._wp
                    do a = 1, nd
                        do b = 1, nd
                            tdotru(a) = tdotru(a) + tb(a, b)*rhou(b)
                        end do
                        q_cls_vf(i_wtu + a)%sf(j, k, l) = real((real(q_stat_vf(lso_stat_rhotau_u_beg + a - 1)%sf(j, k, l), &
                                 & wp) - tdotru(a))/rho_b, stp)
                    end do
                    ! R_mu_sg (filtered minus resolved, same tau convention as the solver)
                    do c = 1, nt
                        tau_res = lso_mu*(dudx(ti(c), tj(c)) + dudx(tj(c), ti(c)))
                        if (ti(c) == tj(c)) tau_res = tau_res - (2._wp/3._wp)*lso_mu*div_u
                        q_cls_vf(i_rmu + c)%sf(j, k, l) = real(tb(ti(c), tj(c)) - tau_res, stp)
                    end do
                    ! R_lam_sg (stored q = -lambda grad T; resolved uses the same sign)
                    do a = 1, nd
                        q_cls_vf(i_rlam + a)%sf(j, k, l) = real(real(q_stat_vf(lso_stat_q_beg + a - 1)%sf(j, k, l), &
                                 & wp)/wloc - (-fluid_k_therm(1)*dT(a)), stp)
                    end do
                    ! T_tilde, u_favre
                    q_cls_vf(i_tt + 1)%sf(j, k, l) = uT_vf(nd + 1)%sf(j, k, l)
                    do a = 1, nd
                        q_cls_vf(i_uf + a)%sf(j, k, l) = uT_vf(a)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        do i = 1, nd + 1
            deallocate (uT_vf(i)%sf)
        end do
        do c = 1, size(q_cls_vf), sys_size
            call s_lso_pp_filter_ghost_refresh(q_cls_vf(c:min(c + sys_size - 1, size(q_cls_vf))), 1)
            if (n > 0) call s_lso_pp_filter_ghost_refresh(q_cls_vf(c:min(c + sys_size - 1, size(q_cls_vf))), 2)
            if (p > 0) call s_lso_pp_filter_ghost_refresh(q_cls_vf(c:min(c + sys_size - 1, size(q_cls_vf))), 3)
        end do

    end subroutine s_compute_lso_closure_fields

    !> Refresh ghosts for direction mpi_dir (1=x, 2=y, 3=z): MPI faces via s_mpi_sendrecv_variables_buffers, BC_GHOST_EXTRAP faces
    !! re-extrapolated from the edge cell, BC_PERIODIC faces wrapped locally.
    impure subroutine s_lso_pp_filter_ghost_refresh(q_cons_vf, mpi_dir)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer, intent(in)               :: mpi_dir
        integer                           :: i, j, k, l, beg_bc, end_bc, nv

        nv = size(q_cons_vf)

        select case (mpi_dir)
        case (1)
            beg_bc = bc_x%beg
            end_bc = bc_x%end
        case (2)
            beg_bc = bc_y%beg
            end_bc = bc_y%end
        case (3)
            beg_bc = bc_z%beg
            end_bc = bc_z%end
        end select

#ifdef MFC_MPI
        if (beg_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, -1, nv)
        if (end_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, 1, nv)
#endif

        select case (mpi_dir)
        case (1)
            do i = 1, nv
                do l = 0, p
                    do k = 0, n
                        do j = 1, lso_pp_filter_radius
                            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(0, k, l)
                            else if (beg_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(m - j + 1, k, l)
                            end if
                            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(m, k, l)
                            else if (end_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(j - 1, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        case (2)
            do i = 1, nv
                do l = 0, p
                    do k = 1, lso_pp_filter_radius
                        do j = 0, m
                            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, 0, l)
                            else if (beg_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, n - k + 1, l)
                            end if
                            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, n, l)
                            else if (end_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, k - 1, l)
                            end if
                        end do
                    end do
                end do
            end do
        case (3)
            do i = 1, nv
                do l = 1, lso_pp_filter_radius
                    do k = 0, n
                        do j = 0, m
                            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, 0)
                            else if (beg_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, p - l + 1)
                            end if
                            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, p)
                            else if (end_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, l - 1)
                            end if
                        end do
                    end do
                end do
            end do
        end select

    end subroutine s_lso_pp_filter_ghost_refresh

    !> Build the LSO stat product fields from the post_process-filtered conserved state. No IB markers in post_process, so phi_p = 0
    !! and gas_mask = 1 everywhere. CPU loops.
    impure subroutine s_compute_lso_pp_stat_fields(q_cons_vf, w_vf)

        type(scalar_field), intent(in) :: q_cons_vf(:)
        !> gas mask (1 fluid, 0 solid): products of the ORIGINAL state are phase-weighted
        type(scalar_field), intent(in), optional :: w_vf(1:1)
        integer                                  :: i, j, k, l
        real(wp)                                 :: g
        real(wp)                                 :: rho, rho_loc
        real(wp)                                 :: mom1, mom2, mom3
        real(wp)                                 :: u1, u2, u3
        real(wp)                                 :: E_loc, ke, e_int, T_loc
        ! Gradient quantities (viscous pass)
        real(wp) :: rho_jm, rho_jp, rho_km, rho_kp, rho_lm, rho_lp
        real(wp) :: u1_jm, u1_jp, u1_km, u1_kp, u1_lm, u1_lp
        real(wp) :: u2_jm, u2_jp, u2_km, u2_kp, u2_lm, u2_lp
        real(wp) :: u3_jm, u3_jp, u3_km, u3_kp, u3_lm, u3_lp
        real(wp) :: T_jm, T_jp, T_km, T_kp, T_lm, T_lp
        real(wp) :: ddx, ddy, ddz
        real(wp) :: du1dx, du1dy, du1dz
        real(wp) :: du2dx, du2dy, du2dz
        real(wp) :: du3dx, du3dy, du3dz
        real(wp) :: dTdx, dTdy, dTdz
        real(wp) :: div_u
        real(wp) :: tau11, tau12, tau13, tau22, tau23, tau33
        real(wp) :: q1, q2, q3

        ! Pass 1: algebraic products. gas_mask = 1 everywhere (no IB).

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp), sgm_eps)
                    mom1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l), wp)
                    if (n > 0) then
                        mom2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l), wp)
                    else
                        mom2 = 0._wp
                    end if
                    if (p > 0) then
                        mom3 = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l), wp)
                    else
                        mom3 = 0._wp
                    end if
                    E_loc = real(q_cons_vf(eqn_idx%E)%sf(j, k, l), wp)

                    u1 = mom1/rho
                    u2 = mom2/rho
                    u3 = mom3/rho
                    ke = 0.5_wp*(mom1**2 + mom2**2 + mom3**2)/rho
                    e_int = (E_loc - ke)/rho
                    T_loc = e_int/(gammas(1)*lso_R_gas)

                    ! phi_p, rho scalar, rhoke scalar, u_p (no IB in post_process filter)
                    q_lso_pp_stat_vf(lso_stat_phi_p_beg)%sf(j, k, l) = 0._stp
                    q_lso_pp_stat_vf(lso_stat_rho_beg)%sf(j, k, l) = real(rho, stp)
                    q_lso_pp_stat_vf(lso_stat_rhoke_beg)%sf(j, k, l) = real((mom1**2 + mom2**2 + mom3**2)/rho, stp)
                    q_lso_pp_stat_vf(lso_stat_up_beg)%sf(j, k, l) = 0._stp
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_up_beg + 1)%sf(j, k, l) = 0._stp
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_up_beg + 2)%sf(j, k, l) = 0._stp

                    ! rho*u
                    q_lso_pp_stat_vf(lso_stat_rhou_beg)%sf(j, k, l) = real(mom1, stp)
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_rhou_beg + 1)%sf(j, k, l) = real(mom2, stp)
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_rhou_beg + 2)%sf(j, k, l) = real(mom3, stp)

                    ! rho*u*u upper triangle
                    q_lso_pp_stat_vf(lso_stat_rhouu_beg)%sf(j, k, l) = real(mom1*u1, stp)
                    if (n > 0) then
                        q_lso_pp_stat_vf(lso_stat_rhouu_beg + 1)%sf(j, k, l) = real(mom1*u2, stp)
                        if (p > 0) then
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 2)%sf(j, k, l) = real(mom1*u3, stp)
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 3)%sf(j, k, l) = real(mom2*u2, stp)
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 4)%sf(j, k, l) = real(mom2*u3, stp)
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 5)%sf(j, k, l) = real(mom3*u3, stp)
                        else
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 2)%sf(j, k, l) = real(mom2*u2, stp)
                        end if
                    end if

                    ! rho*u*|u|^2
                    q_lso_pp_stat_vf(lso_stat_rhouke_beg)%sf(j, k, l) = real(mom1*(u1**2 + u2**2 + u3**2), stp)
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_rhouke_beg + 1)%sf(j, k, l) = real(mom2*(u1**2 + u2**2 + u3**2), stp)
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_rhouke_beg + 2)%sf(j, k, l) = real(mom3*(u1**2 + u2**2 + u3**2), stp)

                    ! rho*u*T
                    q_lso_pp_stat_vf(lso_stat_rhouT_beg)%sf(j, k, l) = real(mom1*T_loc, stp)
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_rhouT_beg + 1)%sf(j, k, l) = real(mom2*T_loc, stp)
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_rhouT_beg + 2)%sf(j, k, l) = real(mom3*T_loc, stp)

                    ! Zero the viscous-pass slots.
                    do i = lso_stat_tau_beg, lso_stat_tau_end
                        q_lso_pp_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                    do i = lso_stat_q_beg, lso_stat_q_end
                        q_lso_pp_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                    do i = lso_stat_rhotau_u_beg, lso_stat_rhotau_u_end
                        q_lso_pp_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                end do
            end do
        end do

        ! Pass 2: viscous stress, heat flux, viscous power flux from centred diffs.
        if (lso_mu > 0._wp .or. fluid_k_therm(1) > 0._wp) then
            if (n == 0) then
                ! 1D: x-gradients only.
                do j = 0, m
                    rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, 0, 0), wp), sgm_eps)

                    rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, 0, 0), wp), sgm_eps)
                    u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, 0, 0), wp)/rho_jm
                    T_jm = (real(q_cons_vf(eqn_idx%E)%sf(j - 1, 0, 0), wp) - 0.5_wp*u1_jm**2*rho_jm)/(rho_jm*gammas(1)*lso_R_gas)

                    rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, 0, 0), wp), sgm_eps)
                    u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, 0, 0), wp)/rho_jp
                    T_jp = (real(q_cons_vf(eqn_idx%E)%sf(j + 1, 0, 0), wp) - 0.5_wp*u1_jp**2*rho_jp)/(rho_jp*gammas(1)*lso_R_gas)

                    ddx = x_cc(j + 1) - x_cc(j - 1)
                    du1dx = (u1_jp - u1_jm)/ddx
                    dTdx = (T_jp - T_jm)/ddx

                    tau11 = lso_mu*(4._wp/3._wp)*du1dx
                    q1 = -fluid_k_therm(1)*dTdx

                    u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, 0, 0), wp)/rho_loc

                    q_lso_pp_stat_vf(lso_stat_tau_beg)%sf(j, 0, 0) = real(tau11, stp)
                    q_lso_pp_stat_vf(lso_stat_q_beg)%sf(j, 0, 0) = real(q1, stp)
                    q_lso_pp_stat_vf(lso_stat_rhotau_u_beg)%sf(j, 0, 0) = real(rho_loc*tau11*u1, stp)
                end do
            else if (p == 0) then
                ! 2D: x- and y-gradients.
                do k = 0, n
                    do j = 0, m
                        rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, 0), wp), sgm_eps)

                        rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, k, 0), wp), sgm_eps)
                        u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, 0), wp)/rho_jm
                        u2_jm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j - 1, k, 0), wp)/rho_jm
                        T_jm = (real(q_cons_vf(eqn_idx%E)%sf(j - 1, k, 0), &
                                & wp) - 0.5_wp*(u1_jm**2 + u2_jm**2)*rho_jm)/(rho_jm*gammas(1)*lso_R_gas)

                        rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, k, 0), wp), sgm_eps)
                        u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, 0), wp)/rho_jp
                        u2_jp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j + 1, k, 0), wp)/rho_jp
                        T_jp = (real(q_cons_vf(eqn_idx%E)%sf(j + 1, k, 0), &
                                & wp) - 0.5_wp*(u1_jp**2 + u2_jp**2)*rho_jp)/(rho_jp*gammas(1)*lso_R_gas)

                        rho_km = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k - 1, 0), wp), sgm_eps)
                        u1_km = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k - 1, 0), wp)/rho_km
                        u2_km = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, 0), wp)/rho_km
                        T_km = (real(q_cons_vf(eqn_idx%E)%sf(j, k - 1, 0), &
                                & wp) - 0.5_wp*(u1_km**2 + u2_km**2)*rho_km)/(rho_km*gammas(1)*lso_R_gas)

                        rho_kp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k + 1, 0), wp), sgm_eps)
                        u1_kp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k + 1, 0), wp)/rho_kp
                        u2_kp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, 0), wp)/rho_kp
                        T_kp = (real(q_cons_vf(eqn_idx%E)%sf(j, k + 1, 0), &
                                & wp) - 0.5_wp*(u1_kp**2 + u2_kp**2)*rho_kp)/(rho_kp*gammas(1)*lso_R_gas)

                        ddx = x_cc(j + 1) - x_cc(j - 1)
                        ddy = y_cc(k + 1) - y_cc(k - 1)

                        du1dx = (u1_jp - u1_jm)/ddx
                        du1dy = (u1_kp - u1_km)/ddy
                        du2dx = (u2_jp - u2_jm)/ddx
                        du2dy = (u2_kp - u2_km)/ddy
                        dTdx = (T_jp - T_jm)/ddx
                        dTdy = (T_kp - T_km)/ddy

                        div_u = du1dx + du2dy
                        tau11 = lso_mu*(2._wp*du1dx - (2._wp/3._wp)*div_u)
                        tau12 = lso_mu*(du1dy + du2dx)
                        tau22 = lso_mu*(2._wp*du2dy - (2._wp/3._wp)*div_u)
                        q1 = -fluid_k_therm(1)*dTdx
                        q2 = -fluid_k_therm(1)*dTdy

                        u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, 0), wp)/rho_loc
                        u2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, 0), wp)/rho_loc

                        q_lso_pp_stat_vf(lso_stat_tau_beg)%sf(j, k, 0) = real(tau11, stp)
                        q_lso_pp_stat_vf(lso_stat_tau_beg + 1)%sf(j, k, 0) = real(tau12, stp)
                        q_lso_pp_stat_vf(lso_stat_tau_beg + 2)%sf(j, k, 0) = real(tau22, stp)
                        q_lso_pp_stat_vf(lso_stat_q_beg)%sf(j, k, 0) = real(q1, stp)
                        q_lso_pp_stat_vf(lso_stat_q_beg + 1)%sf(j, k, 0) = real(q2, stp)
                        q_lso_pp_stat_vf(lso_stat_rhotau_u_beg)%sf(j, k, 0) = real(rho_loc*(tau11*u1 + tau12*u2), stp)
                        q_lso_pp_stat_vf(lso_stat_rhotau_u_beg + 1)%sf(j, k, 0) = real(rho_loc*(tau12*u1 + tau22*u2), stp)
                    end do
                end do
            else
                ! 3D: x-, y-, z-gradients.
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp), sgm_eps)

                            rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, k, l), wp), sgm_eps)
                            u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, l), wp)/rho_jm
                            u2_jm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j - 1, k, l), wp)/rho_jm
                            u3_jm = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j - 1, k, l), wp)/rho_jm
                            T_jm = (real(q_cons_vf(eqn_idx%E)%sf(j - 1, k, l), &
                                    & wp) - 0.5_wp*(u1_jm**2 + u2_jm**2 + u3_jm**2)*rho_jm)/(rho_jm*gammas(1)*lso_R_gas)

                            rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, k, l), wp), sgm_eps)
                            u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, l), wp)/rho_jp
                            u2_jp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j + 1, k, l), wp)/rho_jp
                            u3_jp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j + 1, k, l), wp)/rho_jp
                            T_jp = (real(q_cons_vf(eqn_idx%E)%sf(j + 1, k, l), &
                                    & wp) - 0.5_wp*(u1_jp**2 + u2_jp**2 + u3_jp**2)*rho_jp)/(rho_jp*gammas(1)*lso_R_gas)

                            rho_km = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k - 1, l), wp), sgm_eps)
                            u1_km = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k - 1, l), wp)/rho_km
                            u2_km = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, l), wp)/rho_km
                            u3_km = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k - 1, l), wp)/rho_km
                            T_km = (real(q_cons_vf(eqn_idx%E)%sf(j, k - 1, l), &
                                    & wp) - 0.5_wp*(u1_km**2 + u2_km**2 + u3_km**2)*rho_km)/(rho_km*gammas(1)*lso_R_gas)

                            rho_kp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k + 1, l), wp), sgm_eps)
                            u1_kp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k + 1, l), wp)/rho_kp
                            u2_kp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, l), wp)/rho_kp
                            u3_kp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k + 1, l), wp)/rho_kp
                            T_kp = (real(q_cons_vf(eqn_idx%E)%sf(j, k + 1, l), &
                                    & wp) - 0.5_wp*(u1_kp**2 + u2_kp**2 + u3_kp**2)*rho_kp)/(rho_kp*gammas(1)*lso_R_gas)

                            rho_lm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l - 1), wp), sgm_eps)
                            u1_lm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l - 1), wp)/rho_lm
                            u2_lm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l - 1), wp)/rho_lm
                            u3_lm = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l - 1), wp)/rho_lm
                            T_lm = (real(q_cons_vf(eqn_idx%E)%sf(j, k, l - 1), &
                                    & wp) - 0.5_wp*(u1_lm**2 + u2_lm**2 + u3_lm**2)*rho_lm)/(rho_lm*gammas(1)*lso_R_gas)

                            rho_lp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l + 1), wp), sgm_eps)
                            u1_lp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l + 1), wp)/rho_lp
                            u2_lp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l + 1), wp)/rho_lp
                            u3_lp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l + 1), wp)/rho_lp
                            T_lp = (real(q_cons_vf(eqn_idx%E)%sf(j, k, l + 1), &
                                    & wp) - 0.5_wp*(u1_lp**2 + u2_lp**2 + u3_lp**2)*rho_lp)/(rho_lp*gammas(1)*lso_R_gas)

                            ddx = x_cc(j + 1) - x_cc(j - 1)
                            ddy = y_cc(k + 1) - y_cc(k - 1)
                            ddz = z_cc(l + 1) - z_cc(l - 1)

                            du1dx = (u1_jp - u1_jm)/ddx
                            du1dy = (u1_kp - u1_km)/ddy
                            du1dz = (u1_lp - u1_lm)/ddz
                            du2dx = (u2_jp - u2_jm)/ddx
                            du2dy = (u2_kp - u2_km)/ddy
                            du2dz = (u2_lp - u2_lm)/ddz
                            du3dx = (u3_jp - u3_jm)/ddx
                            du3dy = (u3_kp - u3_km)/ddy
                            du3dz = (u3_lp - u3_lm)/ddz
                            dTdx = (T_jp - T_jm)/ddx
                            dTdy = (T_kp - T_km)/ddy
                            dTdz = (T_lp - T_lm)/ddz

                            div_u = du1dx + du2dy + du3dz
                            tau11 = lso_mu*(2._wp*du1dx - (2._wp/3._wp)*div_u)
                            tau12 = lso_mu*(du1dy + du2dx)
                            tau13 = lso_mu*(du1dz + du3dx)
                            tau22 = lso_mu*(2._wp*du2dy - (2._wp/3._wp)*div_u)
                            tau23 = lso_mu*(du2dz + du3dy)
                            tau33 = lso_mu*(2._wp*du3dz - (2._wp/3._wp)*div_u)
                            q1 = -fluid_k_therm(1)*dTdx
                            q2 = -fluid_k_therm(1)*dTdy
                            q3 = -fluid_k_therm(1)*dTdz

                            u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l), wp)/rho_loc
                            u2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l), wp)/rho_loc
                            u3 = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l), wp)/rho_loc

                            q_lso_pp_stat_vf(lso_stat_tau_beg)%sf(j, k, l) = real(tau11, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 1)%sf(j, k, l) = real(tau12, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 2)%sf(j, k, l) = real(tau13, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 3)%sf(j, k, l) = real(tau22, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 4)%sf(j, k, l) = real(tau23, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 5)%sf(j, k, l) = real(tau33, stp)
                            q_lso_pp_stat_vf(lso_stat_q_beg)%sf(j, k, l) = real(q1, stp)
                            q_lso_pp_stat_vf(lso_stat_q_beg + 1)%sf(j, k, l) = real(q2, stp)
                            q_lso_pp_stat_vf(lso_stat_q_beg + 2)%sf(j, k, l) = real(q3, stp)
                            q_lso_pp_stat_vf(lso_stat_rhotau_u_beg)%sf(j, k, l) = real(rho_loc*(tau11*u1 + tau12*u2 + tau13*u3), &
                                             & stp)
                            q_lso_pp_stat_vf(lso_stat_rhotau_u_beg + 1)%sf(j, k, &
                                             & l) = real(rho_loc*(tau12*u1 + tau22*u2 + tau23*u3), stp)
                            q_lso_pp_stat_vf(lso_stat_rhotau_u_beg + 2)%sf(j, k, &
                                             & l) = real(rho_loc*(tau13*u1 + tau23*u2 + tau33*u3), stp)
                        end do
                    end do
                end do
            end if
        end if

        ! Phase weighting for products formed from the unfiltered state: F[X] = filter(g X), phi_p = filter(1 - g).
        if (present(w_vf)) then
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        g = real(w_vf(1)%sf(j, k, l), wp)
                        q_lso_pp_stat_vf(lso_stat_phi_p_beg)%sf(j, k, l) = real(1._wp - g, stp)
                        do i = lso_stat_phi_p_end + 1, n_lso_stat
                            q_lso_pp_stat_vf(i)%sf(j, k, l) = real(real(q_lso_pp_stat_vf(i)%sf(j, k, l), wp)*g, stp)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_compute_lso_pp_stat_fields

    !> Apply the post_process LSO filter to the stat products in place (chunks of sys_size fields, the size of the shared MPI halo
    !! buffers), through ghost-extended temporaries since q_lso_pp_stat_vf carries no halo.
    impure subroutine s_filter_lso_pp_stat_fields()

        type(scalar_field), allocatable :: q_tmp_vf(:)
        integer                         :: i, c, nc, j, k, l

        allocate (q_tmp_vf(1:sys_size))
        do i = 1, sys_size
            allocate (q_tmp_vf(i)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
        end do
        do c = 1, n_lso_stat, sys_size
            nc = min(sys_size, n_lso_stat - c + 1)
            do i = 1, nc
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_tmp_vf(i)%sf(j, k, l) = q_lso_pp_stat_vf(c + i - 1)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            call s_apply_lso_pp_filter(q_tmp_vf(1:nc))
            do i = 1, nc
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_lso_pp_stat_vf(c + i - 1)%sf(j, k, l) = q_tmp_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        end do
        do i = 1, sys_size
            deallocate (q_tmp_vf(i)%sf)
        end do
        deallocate (q_tmp_vf)

    end subroutine s_filter_lso_pp_stat_fields

end module m_lso_pp_filter
