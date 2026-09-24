!>
!! @file
!! @brief Contains module m_lso_filter

#:include 'macros.fpp'

!> @brief LSO variable-weight Gaussian filter for conserved variables at save steps.
!!
!! 9-point symmetric FIR stencil
!!   H(xi) = a(1) + 2*[a(2)*cos(xi) + a(3)*cos(2*xi) + a(4)*cos(3*xi) + a(5)*cos(4*xi)],
!! composed over several passes to approximate G(xi) = exp(-sigma^2*xi^2/2). Per-pass
!! weights a(1:5) come from the Python BCD design. Applied direction by direction.
module m_lso_filter

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_constants
    use m_eos, only: s_phase_temperature, eoss, isentrope_B
    use m_variables_conversion, only: s_convert_conservative_to_primitive_variables
    use m_ibm, only: ib_markers
    use m_nvtx

    implicit none

    private

    public :: s_initialize_lso_filter_module, s_copy_and_apply_lso_filter, s_lso_stride_sample, s_finalize_lso_filter_module, &
        & q_filt_vf, q_filt_ds_vf, q_lso_stat_vf, q_lso_stat_ds_vf, s_lso_stat_stride_sample, q_lso_mask_vf, q_lso_mask_ds_vf, &
        & s_lso_filter_stage2, s_apply_lso_filter_coarse

    ! Floor on the normalized-convolution denominator filter(m).
    real(wp), parameter :: lso_w_floor = 1.0e-3_wp

    ! Scratch buffer for one directional pass (interior only).
    real(wp), allocatable, dimension(:,:,:) :: lso_tmp
    $:GPU_DECLARE(create='[lso_tmp]')

    ! Filtered copy of the conserved variables (lso_filter_wrt = T).
    type(scalar_field), allocatable :: q_filt_vf(:)

    ! Gas mask (1 in fluid, 0 in solid), the denominator of the normalized convolution filter(m*q)/filter(m).
    type(scalar_field), allocatable :: q_lso_mask_vf(:)

    ! Coarsened version (lso_down_sample_factor > 1).
    type(scalar_field), allocatable :: q_filt_ds_vf(:)

    ! Coarse-grid (stage-2) arrays carry lso_crs_gw ghost layers per active direction; stage 2 runs on the host.
    integer, parameter    :: lso_crs_gw = 4
    integer               :: crs_lo(3), crs_hi(3)
    real(wp), allocatable :: lso2_tmp(:,:,:)
#ifdef MFC_MPI
    real(wp), allocatable :: buff_crs_send(:), buff_crs_recv(:)
#endif

    ! Coarsened filtered gas mask w = filter(m), written so post_process can compose filter2(w*qhat)/filter2(w).
    type(scalar_field), allocatable :: q_lso_mask_ds_vf(:)

    ! Filtered statistical product fields (lso_stat_wrt = T).
    type(scalar_field), allocatable :: q_lso_stat_vf(:)
    type(scalar_field), allocatable :: q_lso_stat_ds_vf(:)
    type(scalar_field), allocatable :: q_lso_T_vf(:)

    ! Host-side halo buffers for the stat exchange (n_lso_stat can exceed sys_size).
#ifdef MFC_MPI
    real(wp), allocatable :: buff_stat_send(:), buff_stat_recv(:)
#endif

contains

    impure subroutine s_initialize_lso_filter_module()

        integer :: i

        @:ALLOCATE(lso_tmp(0:m, 0:n, 0:p))

        if (lso_filter_wrt) then
            @:ALLOCATE(q_filt_vf(1:sys_size))
            do i = 1, sys_size
                @:ALLOCATE(q_filt_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
            do i = 1, sys_size
                @:ACC_SETUP_SFs(q_filt_vf(i))
            end do

            ! Gas-mask denominator field for IB-aware normalized filtering.
            if (ib) then
                @:ALLOCATE(q_lso_mask_vf(1:1))
                @:ALLOCATE(q_lso_mask_vf(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_lso_mask_vf(1))
            end if

            if (lso_down_sample_factor > 1) then
                crs_lo = 0; crs_hi = 0
                crs_lo(1) = -lso_crs_gw; crs_hi(1) = m_lso_ds + lso_crs_gw
                if (n_lso_ds > 0) then
                    crs_lo(2) = -lso_crs_gw; crs_hi(2) = n_lso_ds + lso_crs_gw
                end if
                if (p_lso_ds > 0) then
                    crs_lo(3) = -lso_crs_gw; crs_hi(3) = p_lso_ds + lso_crs_gw
                end if
                @:ALLOCATE(q_filt_ds_vf(1:sys_size))
                do i = 1, sys_size
                    @:ALLOCATE(q_filt_ds_vf(i)%sf(crs_lo(1):crs_hi(1), crs_lo(2):crs_hi(2), crs_lo(3):crs_hi(3)))
                end do
                if (ib) then
                    @:ALLOCATE(q_lso_mask_ds_vf(1:1))
                    @:ALLOCATE(q_lso_mask_ds_vf(1)%sf(crs_lo(1):crs_hi(1), crs_lo(2):crs_hi(2), crs_lo(3):crs_hi(3)))
                end if
                if (lso2_n_passes_x > 0) then
                    allocate (lso2_tmp(0:m_lso_ds,0:n_lso_ds,0:p_lso_ds))
#ifdef MFC_MPI
                    block
                        integer :: nv_max, face_max
                        nv_max = max(sys_size, max(n_lso_stat, 1))
                        face_max = max((n_lso_ds + 1)*(p_lso_ds + 1), (m_lso_ds + 1)*(p_lso_ds + 1), (m_lso_ds + 1)*(n_lso_ds + 1))
                        allocate (buff_crs_send(0:nv_max*lso_crs_gw*face_max - 1))
                        allocate (buff_crs_recv(0:nv_max*lso_crs_gw*face_max - 1))
                    end block
#endif
                end if
            end if

            if (lso_stat_wrt .and. n_lso_stat > 0) then
                @:ALLOCATE(q_lso_T_vf(1:1))
                @:ALLOCATE(q_lso_T_vf(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_lso_T_vf(1))
                @:ALLOCATE(q_lso_stat_vf(1:n_lso_stat))
                do i = 1, n_lso_stat
                    @:ALLOCATE(q_lso_stat_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                end do
                do i = 1, n_lso_stat
                    @:ACC_SETUP_SFs(q_lso_stat_vf(i))
                end do
                if (lso_down_sample_factor > 1) then
                    @:ALLOCATE(q_lso_stat_ds_vf(1:n_lso_stat))
                    do i = 1, n_lso_stat
                        @:ALLOCATE(q_lso_stat_ds_vf(i)%sf(crs_lo(1):crs_hi(1), crs_lo(2):crs_hi(2), crs_lo(3):crs_hi(3)))
                    end do
                end if
#ifdef MFC_MPI
                ! Halo buffers for stat fields - sized for the largest face across all active directions.
                block
                    integer :: halo_stat_size
                    if (n > 0) then
                        if (p > 0) then
                            halo_stat_size = buff_size*n_lso_stat*max((n + 1)*(p + 1), (m + 2*buff_size + 1)*(p + 1), &
                                & (m + 2*buff_size + 1)*(n + 2*buff_size + 1)) - 1
                        else
                            halo_stat_size = buff_size*n_lso_stat*max(n + 1, m + 2*buff_size + 1) - 1
                        end if
                    else
                        halo_stat_size = buff_size*n_lso_stat - 1
                    end if
                    @:ALLOCATE(buff_stat_send(0:halo_stat_size))
                    @:ALLOCATE(buff_stat_recv(0:halo_stat_size))
                end block
#endif
            end if
        end if

    end subroutine s_initialize_lso_filter_module

    impure subroutine s_finalize_lso_filter_module()

        integer :: i

        @:DEALLOCATE(lso_tmp)

        if (lso_filter_wrt) then
            do i = 1, sys_size
                @:DEALLOCATE(q_filt_vf(i)%sf)
            end do
            @:DEALLOCATE(q_filt_vf)

            if (ib) then
                @:DEALLOCATE(q_lso_mask_vf(1)%sf)
                @:DEALLOCATE(q_lso_mask_vf)
            end if

            if (lso_down_sample_factor > 1) then
                do i = 1, sys_size
                    @:DEALLOCATE(q_filt_ds_vf(i)%sf)
                end do
                @:DEALLOCATE(q_filt_ds_vf)
                if (ib) then
                    @:DEALLOCATE(q_lso_mask_ds_vf(1)%sf)
                    @:DEALLOCATE(q_lso_mask_ds_vf)
                end if
                if (allocated(lso2_tmp)) deallocate (lso2_tmp)
#ifdef MFC_MPI
                if (allocated(buff_crs_send)) deallocate (buff_crs_send, buff_crs_recv)
#endif
            end if

            if (lso_stat_wrt .and. n_lso_stat > 0) then
                @:DEALLOCATE(q_lso_T_vf(1)%sf)
                @:DEALLOCATE(q_lso_T_vf)
                if (lso_down_sample_factor > 1) then
                    do i = 1, n_lso_stat
                        @:DEALLOCATE(q_lso_stat_ds_vf(i)%sf)
                    end do
                    @:DEALLOCATE(q_lso_stat_ds_vf)
                end if
                do i = 1, n_lso_stat
                    @:DEALLOCATE(q_lso_stat_vf(i)%sf)
                end do
                @:DEALLOCATE(q_lso_stat_vf)
#ifdef MFC_MPI
                if (allocated(buff_stat_send)) then
                    @:DEALLOCATE(buff_stat_send)
                end if
                if (allocated(buff_stat_recv)) then
                    @:DEALLOCATE(buff_stat_recv)
                end if
#endif
            end if
        end if

    end subroutine s_finalize_lso_filter_module

    !> Copy q_cons_vf into q_filt_vf on the device and filter in place, leaving the original conserved array untouched for the
    !! primary write.
    impure subroutine s_copy_and_apply_lso_filter(q_cons_vf, q_T_sf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        type(scalar_field), intent(in)    :: q_T_sf
        integer                           :: i, j, k, l
        real(wp)                          :: T

        call nvtxStartRange("LSO-FILTER")

        if (lso_stat_wrt .and. n_lso_stat > 0) then
            call nvtxStartRange("LSO-FILTER-STAT")
            ! Reuse the filter workspace for saved-state primitives, never the live RK state.
            if (chemistry) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_lso_T_vf(1)%sf(j, k, l) = q_T_sf%sf(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            call s_convert_conservative_to_primitive_variables(q_cons_vf, q_lso_T_vf(1), q_filt_vf, idwint)
            if (.not. chemistry) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[T]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            call s_lso_temperature(q_cons_vf, q_filt_vf, j, k, l, T)
                            q_lso_T_vf(1)%sf(j, k, l) = real(T, stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if ((viscous .and. lso_mu > 0._wp) .or. fluid_k_therm(1) > 0._wp) then
                do i = 1, num_dims
                    call s_lso_filter_ghost_refresh(q_cons_vf, i)
                    call s_lso_filter_ghost_refresh(q_lso_T_vf, i)
                end do
            end if
            call s_compute_lso_stat_fields(q_cons_vf, q_lso_T_vf(1))
            call s_apply_lso_stat_filter()
#ifndef FRONTIER_UNIFIED
            do i = 1, n_lso_stat
                $:GPU_UPDATE(host='[q_lso_stat_vf(i)%sf]')
            end do
#endif
            call nvtxEndRange
        end if

        call nvtxStartRange("LSO-FILTER-COPY")
        $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_filt_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxEndRange

        if (ib) then
            ! Normalized convolution q_filt = filter(m*q)/filter(m) with the gas mask m (1 in fluid, 0 in solid).
            call nvtxStartRange("LSO-FILTER-MASK")
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ib_markers%sf(j, k, l) > 0) then
                            q_lso_mask_vf(1)%sf(j, k, l) = 0._stp
                        else
                            q_lso_mask_vf(1)%sf(j, k, l) = 1._stp
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_filt_vf(i)%sf(j, k, l) = real(real(q_filt_vf(i)%sf(j, k, l), wp)*real(q_lso_mask_vf(1)%sf(j, k, l), &
                                      & wp), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call nvtxEndRange

            call s_apply_lso_filter(q_filt_vf)
            call s_apply_lso_filter(q_lso_mask_vf)

            ! Normalize everywhere (the solid interior takes the fluid average). Under the two-stage
            ! pyramid the numerator and denominator stay unnormalized until s_lso_filter_stage2.
            if (lso2_n_passes_x <= 0) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_filt_vf(i)%sf(j, k, l) = real(real(q_filt_vf(i)%sf(j, k, l), &
                                          & wp)/max(real(q_lso_mask_vf(1)%sf(j, k, l), wp), lso_w_floor), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else
            call s_apply_lso_filter(q_filt_vf)
        end if

        call nvtxEndRange

    end subroutine s_copy_and_apply_lso_filter

    !> Apply the LSO filter in place to q_cons_vf, direction by direction. Passes loop outside variables so a halo exchange can run
    !! between passes.
    !!
    !! @note With forward Euler (stor == 1) this modifies the live state; with any RK
    !! scheme (stor == 2) only the save copy is touched. IBM runs use RK.
    impure subroutine s_apply_lso_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, ipass, j, k, l, nv
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_cons_vf)

        call nvtxStartRange("LSO-FILTER-X")
        call s_lso_filter_ghost_refresh(q_cons_vf, 1)
        do ipass = 1, lso_n_passes_x
            c0 = lso_a_x(1, ipass)
            c1 = lso_a_x(2, ipass)
            c2 = lso_a_x(3, ipass)
            c3 = lso_a_x(4, ipass)
            c4 = lso_a_x(5, ipass)
            do i = 1, nv
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j - 1, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_cons_vf(i)%sf(j - 2, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_cons_vf(i)%sf(j - 3, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_cons_vf(i)%sf(j - 4, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            if (ipass < lso_n_passes_x) then
                call s_lso_filter_ghost_refresh(q_cons_vf, 1)
            end if
        end do
        call nvtxEndRange

        ! y-direction (2D/3D)
        if (n > 0) then
            call nvtxStartRange("LSO-FILTER-Y")
            call s_lso_filter_ghost_refresh(q_cons_vf, 2)
            do ipass = 1, lso_n_passes_y
                c0 = lso_a_y(1, ipass)
                c1 = lso_a_y(2, ipass)
                c2 = lso_a_y(3, ipass)
                c3 = lso_a_y(4, ipass)
                c4 = lso_a_y(5, ipass)
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k - 1, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_cons_vf(i)%sf(j, k - 2, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 2, l), wp)) + c3*(real(q_cons_vf(i)%sf(j, k - 3, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 3, l), wp)) + c4*(real(q_cons_vf(i)%sf(j, k - 4, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 4, l), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_y) then
                    call s_lso_filter_ghost_refresh(q_cons_vf, 2)
                end if
            end do
            call nvtxEndRange
        end if

        ! z-direction (3D)
        if (p > 0) then
            call nvtxStartRange("LSO-FILTER-Z")
            call s_lso_filter_ghost_refresh(q_cons_vf, 3)
            do ipass = 1, lso_n_passes_z
                c0 = lso_a_z(1, ipass)
                c1 = lso_a_z(2, ipass)
                c2 = lso_a_z(3, ipass)
                c3 = lso_a_z(4, ipass)
                c4 = lso_a_z(5, ipass)
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k, l - 1), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_cons_vf(i)%sf(j, k, l - 2), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 2), wp)) + c3*(real(q_cons_vf(i)%sf(j, k, l - 3), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 3), wp)) + c4*(real(q_cons_vf(i)%sf(j, k, l - 4), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 4), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_z) then
                    call s_lso_filter_ghost_refresh(q_cons_vf, 3)
                end if
            end do
            call nvtxEndRange
        end if

    end subroutine s_apply_lso_filter

    !> Resample q_src_vf into the (pre-allocated) q_dst_vf onto the coarsened grid using trilinear interpolation. Output cell j maps
    !! to source position j*m/m_lso_ds (uniformly spaced from 0 to m), so the first and last cells are always exact and the domain
    !! is fully preserved regardless of divisibility by the stride factor.
    impure subroutine s_lso_stride_sample(q_src_vf, q_dst_vf)

        type(scalar_field), intent(in)    :: q_src_vf(:)
        type(scalar_field), intent(inout) :: q_dst_vf(:)
        integer                           :: i, j, k, l, nv
        integer                           :: j0, j1, k0, k1, l0, l1
        integer                           :: sidx(3)
        real(wp)                          :: alpha, beta, gamma_pos, wj, wk, wl

        ! Global offset of this rank's block (start_idx is only allocated with parallel_io).

        sidx = 0
        if (allocated(start_idx)) sidx(1:size(start_idx)) = start_idx

        nv = size(q_src_vf)
        do i = 1, nv
            do l = 0, p_lso_ds
                ! Sample at the coarse-cell centres: fine position (J + 1/2)(m_glb + 1)/(m_glb_ds + 1) - 1/2 for the
                ! global coarse index J = start_idx/factor + j, so the sample spacing equals the coarse cell size
                ! (interior-only reads for factor >= 2).
                if (p_lso_ds > 0) then
                    gamma_pos = (real(sidx(3)/lso_down_sample_factor + l, wp) + 0.5_wp)*real(p_glb + 1, &
                                 & wp)/real(p_glb_lso_ds + 1, wp) - 0.5_wp - real(sidx(3), wp)
                    l0 = floor(gamma_pos); l1 = l0 + 1; wl = gamma_pos - real(l0, wp)
                else
                    l0 = 0; l1 = 0; wl = 0._wp
                end if

                do k = 0, n_lso_ds
                    if (n_lso_ds > 0) then
                        beta = (real(sidx(2)/lso_down_sample_factor + k, wp) + 0.5_wp)*real(n_glb + 1, wp)/real(n_glb_lso_ds + 1, &
                                & wp) - 0.5_wp - real(sidx(2), wp)
                        k0 = floor(beta); k1 = k0 + 1; wk = beta - real(k0, wp)
                    else
                        k0 = 0; k1 = 0; wk = 0._wp
                    end if

                    do j = 0, m_lso_ds
                        if (m_lso_ds > 0) then
                            alpha = (real(sidx(1)/lso_down_sample_factor + j, wp) + 0.5_wp)*real(m_glb + 1, &
                                     & wp)/real(m_glb_lso_ds + 1, wp) - 0.5_wp - real(sidx(1), wp)
                            j0 = floor(alpha); j1 = j0 + 1; wj = alpha - real(j0, wp)
                        else
                            j0 = 0; j1 = 0; wj = 0._wp
                        end if

                        q_dst_vf(i)%sf(j, k, l) = real((1._wp - wj)*(1._wp - wk)*(1._wp - wl)*real(q_src_vf(i)%sf(j0, k0, l0), &
                                 & wp) + wj*(1._wp - wk)*(1._wp - wl)*real(q_src_vf(i)%sf(j1, k0, l0), &
                                 & wp) + (1._wp - wj)*wk*(1._wp - wl)*real(q_src_vf(i)%sf(j0, k1, l0), &
                                 & wp) + wj*wk*(1._wp - wl)*real(q_src_vf(i)%sf(j1, k1, l0), &
                                 & wp) + (1._wp - wj)*(1._wp - wk)*wl*real(q_src_vf(i)%sf(j0, k0, l1), &
                                 & wp) + wj*(1._wp - wk)*wl*real(q_src_vf(i)%sf(j1, k0, l1), &
                                 & wp) + (1._wp - wj)*wk*wl*real(q_src_vf(i)%sf(j0, k1, l1), &
                                 & wp) + wj*wk*wl*real(q_src_vf(i)%sf(j1, k1, l1), wp), stp)
                    end do
                end do
            end do
        end do

    end subroutine s_lso_stride_sample

    !> Host halo refresh for the coarse (decimated) grid used by the stage-2 pyramid cascade: MPI exchange of lso_crs_gw layers on
    !! decomposed directions (sequential counter packing, identical nesting on pack and unpack), then periodic wrap or edge clamp
    !! for physical boundaries.
    impure subroutine s_lso_coarse_ghost_refresh(q_vf, mpi_dir)

        type(scalar_field), intent(inout) :: q_vf(:)
        integer, intent(in)               :: mpi_dir
        integer                           :: nv, i, j, k, l, g
        integer                           :: beg_bc, end_bc, grid_dim

#ifdef MFC_MPI
        integer :: ierr, r, cnt, pack_offset, unpack_offset, pbc_loc, side
        integer :: dst_proc, src_proc, send_tag, recv_tag
        integer :: beg_end(2)
        logical :: beg_end_geq_0
#endif

        nv = size(q_vf)
        select case (mpi_dir)
        case (1)
            beg_bc = bc_x%beg; end_bc = bc_x%end; grid_dim = m_lso_ds
        case (2)
            beg_bc = bc_y%beg; end_bc = bc_y%end; grid_dim = n_lso_ds
        case (3)
            beg_bc = bc_z%beg; end_bc = bc_z%end; grid_dim = p_lso_ds
        end select

#ifdef MFC_MPI
        beg_end = (/beg_bc, end_bc/)
        select case (mpi_dir)
        case (1)
            cnt = nv*lso_crs_gw*(n_lso_ds + 1)*(p_lso_ds + 1)
        case (2)
            cnt = nv*lso_crs_gw*(m_lso_ds + 1)*(p_lso_ds + 1)
        case (3)
            cnt = nv*lso_crs_gw*(m_lso_ds + 1)*(n_lso_ds + 1)
        end select

        do side = 1, 2
            pbc_loc = 2*side - 3  ! -1 (beg) then +1 (end)
            if (beg_end(side) < 0) cycle
            beg_end_geq_0 = beg_end(max(pbc_loc, 0) - pbc_loc + 1) >= 0
            send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, pbc_loc == 1))
            recv_tag = f_logical_to_int(pbc_loc == 1)
            dst_proc = beg_end(1 + f_logical_to_int(f_xor(pbc_loc == 1, beg_end_geq_0)))
            src_proc = beg_end(1 + f_logical_to_int(pbc_loc == 1))
            pack_offset = 0
            if (f_xor(pbc_loc == 1, beg_end_geq_0)) pack_offset = grid_dim - lso_crs_gw + 1
            unpack_offset = 0
            if (pbc_loc == 1) unpack_offset = grid_dim + lso_crs_gw + 1

            r = -1
            select case (mpi_dir)
            case (1)
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do g = 0, lso_crs_gw - 1
                            do i = 1, nv
                                r = r + 1
                                buff_crs_send(r) = real(q_vf(i)%sf(g + pack_offset, k, l), wp)
                            end do
                        end do
                    end do
                end do
            case (2)
                do l = 0, p_lso_ds
                    do g = 0, lso_crs_gw - 1
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                buff_crs_send(r) = real(q_vf(i)%sf(j, g + pack_offset, l), wp)
                            end do
                        end do
                    end do
                end do
            case (3)
                do g = 0, lso_crs_gw - 1
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                buff_crs_send(r) = real(q_vf(i)%sf(j, k, g + pack_offset), wp)
                            end do
                        end do
                    end do
                end do
            end select

            call MPI_SENDRECV(buff_crs_send, cnt, mpi_p, dst_proc, send_tag, buff_crs_recv, cnt, mpi_p, src_proc, recv_tag, &
                              & MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            r = -1
            select case (mpi_dir)
            case (1)
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do g = -lso_crs_gw, -1
                            do i = 1, nv
                                r = r + 1
                                q_vf(i)%sf(g + unpack_offset, k, l) = real(buff_crs_recv(r), stp)
                            end do
                        end do
                    end do
                end do
            case (2)
                do l = 0, p_lso_ds
                    do g = -lso_crs_gw, -1
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                q_vf(i)%sf(j, g + unpack_offset, l) = real(buff_crs_recv(r), stp)
                            end do
                        end do
                    end do
                end do
            case (3)
                do g = -lso_crs_gw, -1
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                q_vf(i)%sf(j, k, g + unpack_offset) = real(buff_crs_recv(r), stp)
                            end do
                        end do
                    end do
                end do
            end select
        end do
#endif

        ! Physical boundaries: periodic wrap, else clamp to the edge cell.
        do i = 1, nv
            select case (mpi_dir)
            case (1)
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do g = 1, lso_crs_gw
                            if (beg_bc == BC_PERIODIC) then
                                q_vf(i)%sf(-g, k, l) = q_vf(i)%sf(grid_dim - g + 1, k, l)
                            else if (beg_bc < 0) then
                                q_vf(i)%sf(-g, k, l) = q_vf(i)%sf(0, k, l)
                            end if
                            if (end_bc == BC_PERIODIC) then
                                q_vf(i)%sf(grid_dim + g, k, l) = q_vf(i)%sf(g - 1, k, l)
                            else if (end_bc < 0) then
                                q_vf(i)%sf(grid_dim + g, k, l) = q_vf(i)%sf(grid_dim, k, l)
                            end if
                        end do
                    end do
                end do
            case (2)
                do l = 0, p_lso_ds
                    do g = 1, lso_crs_gw
                        do j = 0, m_lso_ds
                            if (beg_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, -g, l) = q_vf(i)%sf(j, grid_dim - g + 1, l)
                            else if (beg_bc < 0) then
                                q_vf(i)%sf(j, -g, l) = q_vf(i)%sf(j, 0, l)
                            end if
                            if (end_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, grid_dim + g, l) = q_vf(i)%sf(j, g - 1, l)
                            else if (end_bc < 0) then
                                q_vf(i)%sf(j, grid_dim + g, l) = q_vf(i)%sf(j, grid_dim, l)
                            end if
                        end do
                    end do
                end do
            case (3)
                do g = 1, lso_crs_gw
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            if (beg_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, k, -g) = q_vf(i)%sf(j, k, grid_dim - g + 1)
                            else if (beg_bc < 0) then
                                q_vf(i)%sf(j, k, -g) = q_vf(i)%sf(j, k, 0)
                            end if
                            if (end_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, k, grid_dim + g) = q_vf(i)%sf(j, k, g - 1)
                            else if (end_bc < 0) then
                                q_vf(i)%sf(j, k, grid_dim + g) = q_vf(i)%sf(j, k, grid_dim)
                            end if
                        end do
                    end do
                end do
            end select
        end do

    end subroutine s_lso_coarse_ghost_refresh

    !> Stage-2 cascade of the two-stage in-situ pyramid: apply the lso2_a_* passes to the coarse-grid fields in place (host).
    impure subroutine s_apply_lso_filter_coarse(q_vf)

        type(scalar_field), intent(inout) :: q_vf(:)
        integer                           :: i, ipass, j, k, l, nv
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_vf)

        call s_lso_coarse_ghost_refresh(q_vf, 1)
        do ipass = 1, lso2_n_passes_x
            c0 = lso2_a_x(1, ipass); c1 = lso2_a_x(2, ipass); c2 = lso2_a_x(3, ipass)
            c3 = lso2_a_x(4, ipass); c4 = lso2_a_x(5, ipass)
            do i = 1, nv
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            lso2_tmp(j, k, l) = c0*real(q_vf(i)%sf(j, k, l), wp) + c1*(real(q_vf(i)%sf(j - 1, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_vf(i)%sf(j - 2, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_vf(i)%sf(j - 3, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_vf(i)%sf(j - 4, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            q_vf(i)%sf(j, k, l) = real(lso2_tmp(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            if (ipass < lso2_n_passes_x) call s_lso_coarse_ghost_refresh(q_vf, 1)
        end do

        if (n_lso_ds > 0) then
            call s_lso_coarse_ghost_refresh(q_vf, 2)
            do ipass = 1, lso2_n_passes_y
                c0 = lso2_a_y(1, ipass); c1 = lso2_a_y(2, ipass); c2 = lso2_a_y(3, ipass)
                c3 = lso2_a_y(4, ipass); c4 = lso2_a_y(5, ipass)
                do i = 1, nv
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                lso2_tmp(j, k, l) = c0*real(q_vf(i)%sf(j, k, l), wp) + c1*(real(q_vf(i)%sf(j, k - 1, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_vf(i)%sf(j, k - 2, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 2, l), wp)) + c3*(real(q_vf(i)%sf(j, k - 3, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 3, l), wp)) + c4*(real(q_vf(i)%sf(j, k - 4, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 4, l), wp))
                            end do
                        end do
                    end do
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                q_vf(i)%sf(j, k, l) = real(lso2_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso2_n_passes_y) call s_lso_coarse_ghost_refresh(q_vf, 2)
            end do
        end if

        if (p_lso_ds > 0) then
            call s_lso_coarse_ghost_refresh(q_vf, 3)
            do ipass = 1, lso2_n_passes_z
                c0 = lso2_a_z(1, ipass); c1 = lso2_a_z(2, ipass); c2 = lso2_a_z(3, ipass)
                c3 = lso2_a_z(4, ipass); c4 = lso2_a_z(5, ipass)
                do i = 1, nv
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                lso2_tmp(j, k, l) = c0*real(q_vf(i)%sf(j, k, l), wp) + c1*(real(q_vf(i)%sf(j, k, l - 1), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_vf(i)%sf(j, k, l - 2), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 2), wp)) + c3*(real(q_vf(i)%sf(j, k, l - 3), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 3), wp)) + c4*(real(q_vf(i)%sf(j, k, l - 4), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 4), wp))
                            end do
                        end do
                    end do
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                q_vf(i)%sf(j, k, l) = real(lso2_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso2_n_passes_z) call s_lso_coarse_ghost_refresh(q_vf, 3)
            end do
        end if

    end subroutine s_apply_lso_filter_coarse

    !> Complete the two-stage in-situ pyramid on the downsampled arrays: stage-2 filter the decimated numerator (and mask
    !! denominator when ib), then normalize. No-op unless lso2 passes are configured.
    impure subroutine s_lso_filter_stage2()

        integer :: i, j, k, l

        if (lso2_n_passes_x <= 0) return

        call s_apply_lso_filter_coarse(q_filt_ds_vf)
        if (ib) then
            call s_apply_lso_filter_coarse(q_lso_mask_ds_vf)
            do i = 1, sys_size
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            q_filt_ds_vf(i)%sf(j, k, l) = real(real(q_filt_ds_vf(i)%sf(j, k, l), &
                                         & wp)/max(real(q_lso_mask_ds_vf(1)%sf(j, k, l), wp), lso_w_floor), stp)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_lso_filter_stage2

    !> Refresh ghost cells between filter passes for direction mpi_dir (1=x, 2=y, 3=z). MPI ghosts come via sendrecv, then
    !! IBM-flagged ghosts are zero-extrapolated to keep particle velocity out of the fluid stencil. BC_GHOST_EXTRAP faces are
    !! re-extrapolated from the current edge cell; other physical BCs are left as is.
    impure subroutine s_lso_filter_ghost_refresh(q_cons_vf, mpi_dir)

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

        ! MPI rank boundaries
#ifdef MFC_MPI
        if (beg_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, -1, nv)
        if (end_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, 1, nv)
#endif

        ! BC_GHOST_EXTRAP: re-extrapolate from the filtered edge cell each pass.
        select case (mpi_dir)
        case (1)
            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (beg_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(m - j + 1, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (end_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(j - 1, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
#ifdef MFC_MPI
#endif
        case (2)
            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (beg_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, n - k + 1, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (end_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, k - 1, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
#ifdef MFC_MPI
#endif
        case (3)
            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (beg_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, p - l + 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (end_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, nv
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, l - 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
#ifdef MFC_MPI
#endif
        end select

    end subroutine s_lso_filter_ghost_refresh

    !> Build the 9/19/31 stat product fields used by the LES closure. Pass 1 fills the algebraic moments (phi_p, u_p, rho*u,
    !! rho*u*u, rho*u*|u|^2, rho*u*T). Pass 2 adds tau_ij, q_i, and (tau*u)_i from centred gradients when either viscosity or
    !! conductivity is active. masked by gas_mask = 1 - phi_p. Both passes run on the GPU.
    subroutine s_lso_temperature(q_cons_vf, q_prim_vf, j, k, l, T)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), intent(in) :: q_cons_vf(:), q_prim_vf(:)
        integer, intent(in)            :: j, k, l
        real(wp), intent(out)          :: T
        real(wp)                       :: alpha_phase, rho_phase, pres

        pres = real(q_prim_vf(eqn_idx%E)%sf(j, k, l), wp)
        if (model_eqns == model_eqns_gamma_law) then
            alpha_phase = 1._wp
        else
            alpha_phase = max(real(q_cons_vf(eqn_idx%E + 1)%sf(j, k, l), wp), sgm_eps)
        end if
        rho_phase = real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp)/alpha_phase
        if (eoss(1) == eos_stiffened_gas .or. eoss(1) == eos_ideal_gas) then
            T = (pres + isentrope_B(1))/(lso_R_gas*rho_phase)
        else
            call s_phase_temperature(rho_phase, pres, 1, T)
        end if

    end subroutine s_lso_temperature

    impure subroutine s_compute_lso_stat_fields(q_cons_vf, q_T_sf)

        type(scalar_field), intent(in) :: q_cons_vf(:)
        type(scalar_field), intent(in) :: q_T_sf
        integer                        :: i, j, k, l, ib_id
        ! Pass 2 loop bounds (narrowed away from physical boundaries, same logic as s_apply_lso_stat_filter).
        integer  :: j_beg_x, j_end_x
        integer  :: k_beg_y, k_end_y
        integer  :: l_beg_z, l_end_z
        real(wp) :: rho, rho_loc
        real(wp) :: mom1, mom2, mom3
        real(wp) :: u1, u2, u3
        real(wp) :: T_loc
        real(wp) :: phi_p, gas_mask
        real(wp) :: up1, up2, up3
        ! Gradient quantities (viscous pass)
        real(wp) :: rho_jm, rho_jp, rho_km, rho_kp, rho_lm, rho_lp
        real(wp) :: u1_jm, u1_jp, u1_km, u1_kp, u1_lm, u1_lp
        real(wp) :: u2_jm, u2_jp, u2_km, u2_kp, u2_lm, u2_lp
        real(wp) :: u3_jm, u3_jp, u3_km, u3_kp, u3_lm, u3_lp
        real(wp) :: T_jm, T_jp, T_km, T_kp, T_lm, T_lp
        real(wp) :: dx, dy, dz
        real(wp) :: du1dx, du1dy, du1dz
        real(wp) :: du2dx, du2dy, du2dz
        real(wp) :: du3dx, du3dy, du3dz
        real(wp) :: dTdx, dTdy, dTdz
        real(wp) :: div_u
        real(wp) :: tau11, tau12, tau13, tau22, tau23, tau33
        real(wp) :: q1, q2, q3

        ! Pass 1: algebraic products.

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, ib_id, phi_p, gas_mask, up1, up2, up3, rho, mom1, mom2, mom3, u1, &
                            & u2, u3, T_loc]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (ib) then
                        if (ib_markers%sf(j, k, l) > 0) then
                            ib_id = ib_markers%sf(j, k, l)
                            phi_p = 1._wp
                            gas_mask = 0._wp
                            up1 = patch_ib(ib_id)%vel(1)
                            up2 = patch_ib(ib_id)%vel(2)
                            up3 = patch_ib(ib_id)%vel(3)
                        else
                            phi_p = 0._wp
                            gas_mask = 1._wp
                            up1 = 0._wp
                            up2 = 0._wp
                            up3 = 0._wp
                        end if
                    else
                        phi_p = 0._wp
                        gas_mask = 1._wp
                        up1 = 0._wp
                        up2 = 0._wp
                        up3 = 0._wp
                    end if

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
                    u1 = mom1/rho
                    u2 = mom2/rho
                    u3 = mom3/rho
                    T_loc = real(q_T_sf%sf(j, k, l), wp)

                    ! phi_p, rho scalar, rhoke scalar, u_p
                    q_lso_stat_vf(lso_stat_phi_p_beg)%sf(j, k, l) = real(phi_p, stp)
                    q_lso_stat_vf(lso_stat_rho_beg)%sf(j, k, l) = real(gas_mask*rho, stp)
                    q_lso_stat_vf(lso_stat_rhoke_beg)%sf(j, k, l) = real(gas_mask*(mom1**2 + mom2**2 + mom3**2)/rho, stp)
                    q_lso_stat_vf(lso_stat_up_beg)%sf(j, k, l) = real(phi_p*up1, stp)
                    if (n > 0) q_lso_stat_vf(lso_stat_up_beg + 1)%sf(j, k, l) = real(phi_p*up2, stp)
                    if (p > 0) q_lso_stat_vf(lso_stat_up_beg + 2)%sf(j, k, l) = real(phi_p*up3, stp)

                    ! rho*u
                    q_lso_stat_vf(lso_stat_rhou_beg)%sf(j, k, l) = real(gas_mask*mom1, stp)
                    if (n > 0) q_lso_stat_vf(lso_stat_rhou_beg + 1)%sf(j, k, l) = real(gas_mask*mom2, stp)
                    if (p > 0) q_lso_stat_vf(lso_stat_rhou_beg + 2)%sf(j, k, l) = real(gas_mask*mom3, stp)

                    ! rho*u*u upper triangle: 1D (11), 2D (11,12,22), 3D (11,12,13,22,23,33)
                    q_lso_stat_vf(lso_stat_rhouu_beg)%sf(j, k, l) = real(gas_mask*mom1*u1, stp)
                    if (n > 0) then
                        q_lso_stat_vf(lso_stat_rhouu_beg + 1)%sf(j, k, l) = real(gas_mask*mom1*u2, stp)
                        if (p > 0) then
                            q_lso_stat_vf(lso_stat_rhouu_beg + 2)%sf(j, k, l) = real(gas_mask*mom1*u3, stp)
                            q_lso_stat_vf(lso_stat_rhouu_beg + 3)%sf(j, k, l) = real(gas_mask*mom2*u2, stp)
                            q_lso_stat_vf(lso_stat_rhouu_beg + 4)%sf(j, k, l) = real(gas_mask*mom2*u3, stp)
                            q_lso_stat_vf(lso_stat_rhouu_beg + 5)%sf(j, k, l) = real(gas_mask*mom3*u3, stp)
                        else
                            q_lso_stat_vf(lso_stat_rhouu_beg + 2)%sf(j, k, l) = real(gas_mask*mom2*u2, stp)
                        end if
                    end if

                    ! rho*u*|u|^2
                    q_lso_stat_vf(lso_stat_rhouke_beg)%sf(j, k, l) = real(gas_mask*mom1*(u1**2 + u2**2 + u3**2), stp)
                    if (n > 0) q_lso_stat_vf(lso_stat_rhouke_beg + 1)%sf(j, k, l) = real(gas_mask*mom2*(u1**2 + u2**2 + u3**2), stp)
                    if (p > 0) q_lso_stat_vf(lso_stat_rhouke_beg + 2)%sf(j, k, l) = real(gas_mask*mom3*(u1**2 + u2**2 + u3**2), stp)

                    ! rho*u*T
                    q_lso_stat_vf(lso_stat_rhouT_beg)%sf(j, k, l) = real(gas_mask*mom1*T_loc, stp)
                    if (n > 0) q_lso_stat_vf(lso_stat_rhouT_beg + 1)%sf(j, k, l) = real(gas_mask*mom2*T_loc, stp)
                    if (p > 0) q_lso_stat_vf(lso_stat_rhouT_beg + 2)%sf(j, k, l) = real(gas_mask*mom3*T_loc, stp)

                    ! Zero the viscous-pass slots.
                    do i = lso_stat_tau_beg, lso_stat_tau_end
                        q_lso_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                    do i = lso_stat_q_beg, lso_stat_q_end
                        q_lso_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                    do i = lso_stat_rhotau_u_beg, lso_stat_rhotau_u_end
                        q_lso_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Pass 2: viscous stress, heat flux, viscous power flux from centred diffs, phase-weighted by the
        ! gas mask like the algebraic blocks (post_process recovers the intensive averages by dividing by
        ! the filtered mask w). Stencils read the raw (unmasked) state, so the image-point fill inside the
        ! IB supplies the wall condition for surface-adjacent gradients.
        ! Boundary ranges: skip cells within buff_size of a physical boundary, since the centred-difference
        ! stencil would reach into solver BC ghost cells (e.g. slip-wall reflections, inlet values) that are
        ! not appropriate for viscous gradient evaluation. Cells outside these ranges keep tau=q=0 from Pass 1.
        j_beg_x = 0; j_end_x = m
        if (bc_x%beg < 0 .and. bc_x%beg /= BC_PERIODIC) j_beg_x = buff_size
        if (bc_x%end < 0 .and. bc_x%end /= BC_PERIODIC) j_end_x = m - buff_size
        k_beg_y = 0; k_end_y = n
        if (bc_y%beg < 0 .and. bc_y%beg /= BC_PERIODIC) k_beg_y = buff_size
        if (bc_y%end < 0 .and. bc_y%end /= BC_PERIODIC) k_end_y = n - buff_size
        l_beg_z = 0; l_end_z = p
        if (bc_z%beg < 0 .and. bc_z%beg /= BC_PERIODIC) l_beg_z = buff_size
        if (bc_z%end < 0 .and. bc_z%end /= BC_PERIODIC) l_end_z = p - buff_size

        if ((viscous .and. lso_mu > 0._wp) .or. fluid_k_therm(1) > 0._wp) then
            if (n == 0) then
                ! 1D: x-gradients only.
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, rho_loc, gas_mask, rho_jm, rho_jp, u1_jm, u1_jp, T_jm, T_jp, &
                                    & dx, du1dx, dTdx, tau11, q1, u1]')
                do l = 0, p
                    do k = 0, n
                        do j = j_beg_x, j_end_x
                            rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, 0, 0), wp), sgm_eps)
                            gas_mask = 1._wp - real(q_lso_stat_vf(lso_stat_phi_p_beg)%sf(j, 0, 0), wp)

                            rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, 0, 0), wp), sgm_eps)
                            u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, 0, 0), wp)/rho_jm
                            T_jm = real(q_T_sf%sf(j - 1, 0, 0), wp)

                            rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, 0, 0), wp), sgm_eps)
                            u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, 0, 0), wp)/rho_jp
                            T_jp = real(q_T_sf%sf(j + 1, 0, 0), wp)

                            dx = x_cc(j + 1) - x_cc(j - 1)
                            du1dx = (u1_jp - u1_jm)/dx
                            dTdx = (T_jp - T_jm)/dx

                            tau11 = lso_mu*(4._wp/3._wp)*du1dx
                            q1 = -fluid_k_therm(1)*dTdx

                            u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, 0, 0), wp)/rho_loc

                            q_lso_stat_vf(lso_stat_tau_beg)%sf(j, 0, 0) = real(gas_mask*tau11, stp)
                            q_lso_stat_vf(lso_stat_q_beg)%sf(j, 0, 0) = real(gas_mask*q1, stp)
                            q_lso_stat_vf(lso_stat_rhotau_u_beg)%sf(j, 0, 0) = real(gas_mask*rho_loc*tau11*u1, stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (p == 0) then
                ! 2D: x- and y-gradients.
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, rho_loc, gas_mask, rho_jm, rho_jp, rho_km, rho_kp, u1_jm, &
                                    & u1_jp, u1_km, u1_kp, u2_jm, u2_jp, u2_km, u2_kp, T_jm, T_jp, T_km, T_kp, dx, dy, du1dx, &
                                    & du1dy, du2dx, du2dy, dTdx, dTdy, div_u, tau11, tau12, tau22, q1, q2, u1, u2]')
                do l = 0, p
                    do k = k_beg_y, k_end_y
                        do j = j_beg_x, j_end_x
                            rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, 0), wp), sgm_eps)
                            gas_mask = 1._wp - real(q_lso_stat_vf(lso_stat_phi_p_beg)%sf(j, k, 0), wp)

                            rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, k, 0), wp), sgm_eps)
                            u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, 0), wp)/rho_jm
                            u2_jm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j - 1, k, 0), wp)/rho_jm
                            T_jm = real(q_T_sf%sf(j - 1, k, 0), wp)

                            rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, k, 0), wp), sgm_eps)
                            u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, 0), wp)/rho_jp
                            u2_jp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j + 1, k, 0), wp)/rho_jp
                            T_jp = real(q_T_sf%sf(j + 1, k, 0), wp)

                            rho_km = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k - 1, 0), wp), sgm_eps)
                            u1_km = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k - 1, 0), wp)/rho_km
                            u2_km = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, 0), wp)/rho_km
                            T_km = real(q_T_sf%sf(j, k - 1, 0), wp)

                            rho_kp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k + 1, 0), wp), sgm_eps)
                            u1_kp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k + 1, 0), wp)/rho_kp
                            u2_kp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, 0), wp)/rho_kp
                            T_kp = real(q_T_sf%sf(j, k + 1, 0), wp)

                            dx = x_cc(j + 1) - x_cc(j - 1)
                            dy = y_cc(k + 1) - y_cc(k - 1)

                            du1dx = (u1_jp - u1_jm)/dx
                            du1dy = (u1_kp - u1_km)/dy
                            du2dx = (u2_jp - u2_jm)/dx
                            du2dy = (u2_kp - u2_km)/dy
                            dTdx = (T_jp - T_jm)/dx
                            dTdy = (T_kp - T_km)/dy

                            div_u = du1dx + du2dy
                            tau11 = lso_mu*(2._wp*du1dx - (2._wp/3._wp)*div_u)
                            tau12 = lso_mu*(du1dy + du2dx)
                            tau22 = lso_mu*(2._wp*du2dy - (2._wp/3._wp)*div_u)
                            q1 = -fluid_k_therm(1)*dTdx
                            q2 = -fluid_k_therm(1)*dTdy

                            u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, 0), wp)/rho_loc
                            u2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, 0), wp)/rho_loc

                            q_lso_stat_vf(lso_stat_tau_beg)%sf(j, k, 0) = real(gas_mask*tau11, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 1)%sf(j, k, 0) = real(gas_mask*tau12, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 2)%sf(j, k, 0) = real(gas_mask*tau22, stp)
                            q_lso_stat_vf(lso_stat_q_beg)%sf(j, k, 0) = real(gas_mask*q1, stp)
                            q_lso_stat_vf(lso_stat_q_beg + 1)%sf(j, k, 0) = real(gas_mask*q2, stp)
                            q_lso_stat_vf(lso_stat_rhotau_u_beg)%sf(j, k, 0) = real(gas_mask*rho_loc*(tau11*u1 + tau12*u2), stp)
                            q_lso_stat_vf(lso_stat_rhotau_u_beg + 1)%sf(j, k, 0) = real(gas_mask*rho_loc*(tau12*u1 + tau22*u2), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                ! 3D: x-, y-, and z-gradients.
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, rho_loc, gas_mask, rho_jm, rho_jp, rho_km, rho_kp, rho_lm, &
                                    & rho_lp, u1_jm, u1_jp, u1_km, u1_kp, u1_lm, u1_lp, u2_jm, u2_jp, u2_km, u2_kp, u2_lm, u2_lp, &
                                    & u3_jm, u3_jp, u3_km, u3_kp, u3_lm, u3_lp, T_jm, T_jp, T_km, T_kp, T_lm, T_lp, dx, dy, dz, &
                                    & du1dx, du1dy, du1dz, du2dx, du2dy, du2dz, du3dx, du3dy, du3dz, dTdx, dTdy, dTdz, div_u, &
                                    & tau11, tau12, tau13, tau22, tau23, tau33, q1, q2, q3, u1, u2, u3]')
                do l = l_beg_z, l_end_z
                    do k = k_beg_y, k_end_y
                        do j = j_beg_x, j_end_x
                            rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp), sgm_eps)
                            gas_mask = 1._wp - real(q_lso_stat_vf(lso_stat_phi_p_beg)%sf(j, k, l), wp)

                            rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, k, l), wp), sgm_eps)
                            u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, l), wp)/rho_jm
                            u2_jm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j - 1, k, l), wp)/rho_jm
                            u3_jm = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j - 1, k, l), wp)/rho_jm
                            T_jm = real(q_T_sf%sf(j - 1, k, l), wp)

                            rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, k, l), wp), sgm_eps)
                            u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, l), wp)/rho_jp
                            u2_jp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j + 1, k, l), wp)/rho_jp
                            u3_jp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j + 1, k, l), wp)/rho_jp
                            T_jp = real(q_T_sf%sf(j + 1, k, l), wp)

                            rho_km = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k - 1, l), wp), sgm_eps)
                            u1_km = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k - 1, l), wp)/rho_km
                            u2_km = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, l), wp)/rho_km
                            u3_km = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k - 1, l), wp)/rho_km
                            T_km = real(q_T_sf%sf(j, k - 1, l), wp)

                            rho_kp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k + 1, l), wp), sgm_eps)
                            u1_kp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k + 1, l), wp)/rho_kp
                            u2_kp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, l), wp)/rho_kp
                            u3_kp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k + 1, l), wp)/rho_kp
                            T_kp = real(q_T_sf%sf(j, k + 1, l), wp)

                            rho_lm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l - 1), wp), sgm_eps)
                            u1_lm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l - 1), wp)/rho_lm
                            u2_lm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l - 1), wp)/rho_lm
                            u3_lm = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l - 1), wp)/rho_lm
                            T_lm = real(q_T_sf%sf(j, k, l - 1), wp)

                            rho_lp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l + 1), wp), sgm_eps)
                            u1_lp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l + 1), wp)/rho_lp
                            u2_lp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l + 1), wp)/rho_lp
                            u3_lp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l + 1), wp)/rho_lp
                            T_lp = real(q_T_sf%sf(j, k, l + 1), wp)

                            dx = x_cc(j + 1) - x_cc(j - 1)
                            dy = y_cc(k + 1) - y_cc(k - 1)
                            dz = z_cc(l + 1) - z_cc(l - 1)

                            du1dx = (u1_jp - u1_jm)/dx
                            du1dy = (u1_kp - u1_km)/dy
                            du1dz = (u1_lp - u1_lm)/dz
                            du2dx = (u2_jp - u2_jm)/dx
                            du2dy = (u2_kp - u2_km)/dy
                            du2dz = (u2_lp - u2_lm)/dz
                            du3dx = (u3_jp - u3_jm)/dx
                            du3dy = (u3_kp - u3_km)/dy
                            du3dz = (u3_lp - u3_lm)/dz
                            dTdx = (T_jp - T_jm)/dx
                            dTdy = (T_kp - T_km)/dy
                            dTdz = (T_lp - T_lm)/dz

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

                            q_lso_stat_vf(lso_stat_tau_beg)%sf(j, k, l) = real(gas_mask*tau11, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 1)%sf(j, k, l) = real(gas_mask*tau12, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 2)%sf(j, k, l) = real(gas_mask*tau13, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 3)%sf(j, k, l) = real(gas_mask*tau22, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 4)%sf(j, k, l) = real(gas_mask*tau23, stp)
                            q_lso_stat_vf(lso_stat_tau_beg + 5)%sf(j, k, l) = real(gas_mask*tau33, stp)
                            q_lso_stat_vf(lso_stat_q_beg)%sf(j, k, l) = real(gas_mask*q1, stp)
                            q_lso_stat_vf(lso_stat_q_beg + 1)%sf(j, k, l) = real(gas_mask*q2, stp)
                            q_lso_stat_vf(lso_stat_q_beg + 2)%sf(j, k, l) = real(gas_mask*q3, stp)
                            ! (tau u)_i = sum_j tau_ij u_j
                            q_lso_stat_vf(lso_stat_rhotau_u_beg)%sf(j, k, &
                                          & l) = real(gas_mask*rho_loc*(tau11*u1 + tau12*u2 + tau13*u3), stp)
                            q_lso_stat_vf(lso_stat_rhotau_u_beg + 1)%sf(j, k, &
                                          & l) = real(gas_mask*rho_loc*(tau12*u1 + tau22*u2 + tau23*u3), stp)
                            q_lso_stat_vf(lso_stat_rhotau_u_beg + 2)%sf(j, k, &
                                          & l) = real(gas_mask*rho_loc*(tau13*u1 + tau23*u2 + tau33*u3), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_compute_lso_stat_fields

    !> LSO filter pass on q_lso_stat_vf (analogous to s_apply_lso_filter for q_cons_vf). Stencil runs on the GPU; the inter-pass
    !! halo exchange goes through the host because n_lso_stat does not fit the global MPI buffers. Use the same interior extent and
    !! physical-boundary extrapolation as the flow and mask filters to preserve complementary volume fractions.
    impure subroutine s_apply_lso_stat_filter()

        integer  :: i, ipass, j, k, l
        real(wp) :: c0, c1, c2, c3, c4

        call s_lso_stat_ghost_refresh(1)
        do ipass = 1, lso_n_passes_x
            c0 = lso_a_x(1, ipass)
            c1 = lso_a_x(2, ipass)
            c2 = lso_a_x(3, ipass)
            c3 = lso_a_x(4, ipass)
            c4 = lso_a_x(5, ipass)
            do i = 1, n_lso_stat
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_tmp(j, k, l) = c0*real(q_lso_stat_vf(i)%sf(j, k, l), wp) + c1*(real(q_lso_stat_vf(i)%sf(j - 1, k, &
                                    & l), wp) + real(q_lso_stat_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_lso_stat_vf(i)%sf(j - 2, &
                                    & k, l), wp) + real(q_lso_stat_vf(i)%sf(j + 2, k, l), &
                                    & wp)) + c3*(real(q_lso_stat_vf(i)%sf(j - 3, k, l), wp) + real(q_lso_stat_vf(i)%sf(j + 3, k, &
                                    & l), wp)) + c4*(real(q_lso_stat_vf(i)%sf(j - 4, k, l), wp) + real(q_lso_stat_vf(i)%sf(j + 4, &
                                    & k, l), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_lso_stat_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            if (ipass < lso_n_passes_x) then
                call s_lso_stat_ghost_refresh(1)
            end if
        end do

        if (n > 0) then
            call s_lso_stat_ghost_refresh(2)
            do ipass = 1, lso_n_passes_y
                c0 = lso_a_y(1, ipass)
                c1 = lso_a_y(2, ipass)
                c2 = lso_a_y(3, ipass)
                c3 = lso_a_y(4, ipass)
                c4 = lso_a_y(5, ipass)
                do i = 1, n_lso_stat
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_lso_stat_vf(i)%sf(j, k, l), wp) + c1*(real(q_lso_stat_vf(i)%sf(j, &
                                        & k - 1, l), wp) + real(q_lso_stat_vf(i)%sf(j, k + 1, l), &
                                        & wp)) + c2*(real(q_lso_stat_vf(i)%sf(j, k - 2, l), wp) + real(q_lso_stat_vf(i)%sf(j, &
                                        & k + 2, l), wp)) + c3*(real(q_lso_stat_vf(i)%sf(j, k - 3, l), &
                                        & wp) + real(q_lso_stat_vf(i)%sf(j, k + 3, l), wp)) + c4*(real(q_lso_stat_vf(i)%sf(j, &
                                        & k - 4, l), wp) + real(q_lso_stat_vf(i)%sf(j, k + 4, l), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_y) then
                    call s_lso_stat_ghost_refresh(2)
                end if
            end do
        end if

        if (p > 0) then
            call s_lso_stat_ghost_refresh(3)
            do ipass = 1, lso_n_passes_z
                c0 = lso_a_z(1, ipass)
                c1 = lso_a_z(2, ipass)
                c2 = lso_a_z(3, ipass)
                c3 = lso_a_z(4, ipass)
                c4 = lso_a_z(5, ipass)
                do i = 1, n_lso_stat
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_lso_stat_vf(i)%sf(j, k, l), wp) + c1*(real(q_lso_stat_vf(i)%sf(j, k, &
                                        & l - 1), wp) + real(q_lso_stat_vf(i)%sf(j, k, l + 1), &
                                        & wp)) + c2*(real(q_lso_stat_vf(i)%sf(j, k, l - 2), wp) + real(q_lso_stat_vf(i)%sf(j, k, &
                                        & l + 2), wp)) + c3*(real(q_lso_stat_vf(i)%sf(j, k, l - 3), &
                                        & wp) + real(q_lso_stat_vf(i)%sf(j, k, l + 3), wp)) + c4*(real(q_lso_stat_vf(i)%sf(j, k, &
                                        & l - 4), wp) + real(q_lso_stat_vf(i)%sf(j, k, l + 4), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_z) then
                    call s_lso_stat_ghost_refresh(3)
                end if
            end do
        end if

    end subroutine s_apply_lso_stat_filter

    !> Ghost-cell refresh for q_lso_stat_vf in direction mpi_dir (1=x, 2=y, 3=z). Halos are packed on the GPU into buff_stat_*, only
    !! the small halo buffers are transferred host<->device for MPI, BC_GHOST_EXTRAP / BC_PERIODIC are reapplied on the GPU, and the
    !! stat fields never leave device memory.
    impure subroutine s_lso_stat_ghost_refresh(mpi_dir)

        integer, intent(in) :: mpi_dir
        integer             :: i, j, k, l, beg_bc, end_bc

#ifdef MFC_MPI
        integer :: ierr, r, buffer_count, pack_offset, unpack_offset
        integer :: dst_proc, src_proc, send_tag, recv_tag
        integer :: grid_dim
        integer :: beg_end(2)
        logical :: beg_end_geq_0
#endif

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
        ! Dedicated halo exchange (buffers sized for n_lso_stat, packed on the host); tag/proc logic
        ! mirrors s_mpi_sendrecv_variables_buffers.
        select case (mpi_dir)
        case (1)
            beg_end = (/bc_x%beg, bc_x%end/)
            buffer_count = buff_size*n_lso_stat*(n + 1)*(p + 1)
            grid_dim = m
        case (2)
            beg_end = (/bc_y%beg, bc_y%end/)
            buffer_count = buff_size*n_lso_stat*(m + 2*buff_size + 1)*(p + 1)
            grid_dim = n
        case (3)
            beg_end = (/bc_z%beg, bc_z%end/)
            buffer_count = buff_size*n_lso_stat*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)
            grid_dim = p
        end select

        ! pbc_loc = -1 then +1.
        if (beg_bc >= 0) then
            ! pbc_loc = -1
            beg_end_geq_0 = beg_end(2) >= 0
            send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, .false.))
            recv_tag = 0
            dst_proc = beg_end(1 + f_logical_to_int(f_xor(.false., beg_end_geq_0)))
            src_proc = beg_end(1)
            pack_offset = 0
            if (beg_end_geq_0) pack_offset = grid_dim - buff_size + 1
            unpack_offset = 0

            select case (mpi_dir)
            case (1)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, buff_size - 1
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*(j + buff_size*(k + (n + 1)*l))
                                buff_stat_send(r) = real(q_lso_stat_vf(i)%sf(j + pack_offset, k, l), wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (2)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = 0, buff_size - 1
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*(k + buff_size*l))
                                buff_stat_send(r) = real(q_lso_stat_vf(i)%sf(j, k + pack_offset, l), wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (3)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, buff_size - 1
                    do k = -buff_size, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*((k + buff_size) + (n &
                                     & + 2*buff_size + 1)*l))
                                buff_stat_send(r) = real(q_lso_stat_vf(i)%sf(j, k, l + pack_offset), wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end select

#ifndef FRONTIER_UNIFIED
            $:GPU_UPDATE(host='[buff_stat_send]')
#endif
            call MPI_SENDRECV(buff_stat_send, buffer_count, mpi_p, dst_proc, send_tag, buff_stat_recv, buffer_count, mpi_p, &
                              & src_proc, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#ifndef FRONTIER_UNIFIED
            $:GPU_UPDATE(device='[buff_stat_recv]')
#endif

            select case (mpi_dir)
            case (1)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size, -1
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*(j + buff_size*((k + 1) + (n + 1)*l))
                                q_lso_stat_vf(i)%sf(j + unpack_offset, k, l) = real(buff_stat_recv(r), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (2)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = -buff_size, -1
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*(k + buff_size*(1 + l)))
                                q_lso_stat_vf(i)%sf(j, k + unpack_offset, l) = real(buff_stat_recv(r), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (3)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = -buff_size, -1
                    do k = -buff_size, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*((k + buff_size) + (n &
                                     & + 2*buff_size + 1)*(l + buff_size)))
                                q_lso_stat_vf(i)%sf(j, k, l + unpack_offset) = real(buff_stat_recv(r), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end select
        end if

        if (end_bc >= 0) then
            ! pbc_loc = +1
            beg_end_geq_0 = beg_end(1) >= 0
            send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, .true.))
            recv_tag = 1
            dst_proc = beg_end(1 + f_logical_to_int(f_xor(.true., beg_end_geq_0)))
            src_proc = beg_end(2)
            pack_offset = 0
            if (f_xor(.true., beg_end_geq_0)) pack_offset = grid_dim - buff_size + 1
            unpack_offset = grid_dim + buff_size + 1

            select case (mpi_dir)
            case (1)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, buff_size - 1
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*(j + buff_size*(k + (n + 1)*l))
                                buff_stat_send(r) = real(q_lso_stat_vf(i)%sf(j + pack_offset, k, l), wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (2)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = 0, buff_size - 1
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*(k + buff_size*l))
                                buff_stat_send(r) = real(q_lso_stat_vf(i)%sf(j, k + pack_offset, l), wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (3)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, buff_size - 1
                    do k = -buff_size, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*((k + buff_size) + (n &
                                     & + 2*buff_size + 1)*l))
                                buff_stat_send(r) = real(q_lso_stat_vf(i)%sf(j, k, l + pack_offset), wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end select

#ifndef FRONTIER_UNIFIED
            $:GPU_UPDATE(host='[buff_stat_send]')
#endif
            call MPI_SENDRECV(buff_stat_send, buffer_count, mpi_p, dst_proc, send_tag, buff_stat_recv, buffer_count, mpi_p, &
                              & src_proc, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#ifndef FRONTIER_UNIFIED
            $:GPU_UPDATE(device='[buff_stat_recv]')
#endif

            select case (mpi_dir)
            case (1)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size, -1
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*(j + buff_size*((k + 1) + (n + 1)*l))
                                q_lso_stat_vf(i)%sf(j + unpack_offset, k, l) = real(buff_stat_recv(r), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (2)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = 0, p
                    do k = -buff_size, -1
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*(k + buff_size*(1 + l)))
                                q_lso_stat_vf(i)%sf(j, k + unpack_offset, l) = real(buff_stat_recv(r), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            case (3)
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, r]')
                do l = -buff_size, -1
                    do k = -buff_size, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, n_lso_stat
                                r = (i - 1) + n_lso_stat*((j + buff_size) + (m + 2*buff_size + 1)*((k + buff_size) + (n &
                                     & + 2*buff_size + 1)*(l + buff_size)))
                                q_lso_stat_vf(i)%sf(j, k, l + unpack_offset) = real(buff_stat_recv(r), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end select
        end if
#endif

        ! BC_GHOST_EXTRAP: copy the edge cell into the ghost layer.
        ! BC_PERIODIC (single-rank, no MPI partner): wrap the local domain.
        select case (mpi_dir)
        case (1)
            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_lso_stat_vf(i)%sf(-j, k, l) = q_lso_stat_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (beg_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_lso_stat_vf(i)%sf(-j, k, l) = q_lso_stat_vf(i)%sf(m - j + 1, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_lso_stat_vf(i)%sf(m + j, k, l) = q_lso_stat_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (end_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_lso_stat_vf(i)%sf(m + j, k, l) = q_lso_stat_vf(i)%sf(j - 1, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        case (2)
            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, -k, l) = q_lso_stat_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (beg_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, -k, l) = q_lso_stat_vf(i)%sf(j, n - k + 1, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, n + k, l) = q_lso_stat_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (end_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, n + k, l) = q_lso_stat_vf(i)%sf(j, k - 1, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        case (3)
            if (beg_bc == BC_GHOST_EXTRAP .or. (beg_bc < 0 .and. beg_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, k, -l) = q_lso_stat_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (beg_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, k, -l) = q_lso_stat_vf(i)%sf(j, k, p - l + 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (end_bc == BC_GHOST_EXTRAP .or. (end_bc < 0 .and. end_bc /= BC_PERIODIC)) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, k, p + l) = q_lso_stat_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (end_bc == BC_PERIODIC) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = 1, n_lso_stat
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_lso_stat_vf(i)%sf(j, k, p + l) = q_lso_stat_vf(i)%sf(j, k, l - 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end select

    end subroutine s_lso_stat_ghost_refresh

    !> Decimate the stat fields onto the coarse grid (same mapping as s_lso_stride_sample).
    impure subroutine s_lso_stat_stride_sample()

        call s_lso_stride_sample(q_lso_stat_vf(1:n_lso_stat), q_lso_stat_ds_vf(1:n_lso_stat))

    end subroutine s_lso_stat_stride_sample

end module m_lso_filter
