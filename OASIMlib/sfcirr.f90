submodule (oasim) oasim_sfcirr
    use oasim_device, only: light, clrtrans, slingo
#ifdef _OPENACC
    use openacc 
    use nvtx
#endif
    USE OMP_LIB
    implicit none
contains
    module subroutine sfcirr(self, iday, sec_c, slp, wsm, oz, wv, rh, &
                            taua, asymp, ssalb, ccov, rlwp, cdre, error)
        implicit none
        class(calc_unit) :: self
        integer, intent(in) :: iday
        real(kind=real_kind), intent(in) :: sec_c
        real(kind=real_kind), dimension(:), intent(in) :: slp, wsm, oz, wv, rh
        real(kind=real_kind), dimension(:,:), intent(in) :: taua, asymp, ssalb
        real(kind=real_kind), dimension(:), intent(in) :: ccov, rlwp, cdre
        logical, intent(out) :: error

        real(kind=real_kind), parameter :: daypersec = 1.0d0 / 86400.0d0
        real(kind=real_kind) :: rday, daycor, sunz, cosunz, pres, ws, ozone, wvapor, relhum
        real(kind=real_kind) :: cov, clwp, re
        real(kind=real_kind) :: am, vi
        integer :: j, rows_size

        ! Local copies of frequently accessed derived type components
        real(kind=real_kind), dimension(:), allocatable :: fobar, oza, awv, ao, aco2, tab2
        real(kind=real_kind), dimension(:), allocatable :: asl, bsl, csl, dsl, esl, fsl
        
        integer, allocatable :: daylight_idx(:)
        integer :: nvalid, i 
        integer :: ngpus, g, devtype, start_idx, end_idx, chunk
        integer :: local_chunk_size
        
        ! Per-GPU allocatable arrays
        real(kind=real_kind), dimension(:), allocatable :: ed_local, es_local
        real(kind=real_kind), dimension(:), allocatable :: ta_local, wa_local, asym_local
        real(kind=real_kind), dimension(:), allocatable :: rlamu_local, td_local, ts_local
        real(kind=real_kind), dimension(:), allocatable :: tcd_local, tcs_local, tgas_local
        real(kind=real_kind), dimension(:), allocatable :: edclr_local, esclr_local
        real(kind=real_kind), dimension(:), allocatable :: edcld_local, escld_local
        real(kind=real_kind), dimension(:, :), allocatable :: ed_host, es_host

        am = self%lib%init_parameters%am
        vi = self%lib%init_parameters%vi
        rows_size = self%lib%rows
        
        ! Initialize output arrays
        self%eda = 0.0d0
        self%esa = 0.0d0
        
        rday = real(iday, real_kind) + sec_c * daypersec
        daycor = 1.0 + 1.67d-2 * cos(pi2 * (rday - 3.0d0) / 365.0d0)
        daycor = daycor * daycor

        ! Allocate local arrays
        allocate(fobar(rows_size), oza(rows_size), awv(rows_size), ao(rows_size))
        allocate(aco2(rows_size), tab2(rows_size))
        allocate(asl(rows_size), bsl(rows_size), csl(rows_size), dsl(rows_size))
        allocate(esl(rows_size), fsl(rows_size))
        allocate(ed_host(self%p_size, rows_size))
        allocate(es_host(self%p_size, rows_size))

        ! Copy data from derived types to local arrays
        fobar = self%lib%atmo_adapted%tab(:,1)
        oza = self%lib%atmo_adapted%tab(:,3)
        awv = self%lib%atmo_adapted%tab(:,4)
        ao = self%lib%atmo_adapted%tab(:,5)
        aco2 = self%lib%atmo_adapted%tab(:,6)
        tab2 = self%lib%atmo_adapted%tab(:,2)
        
        asl = self%lib%slingo_adapted%tab(:,1)
        bsl = self%lib%slingo_adapted%tab(:,2)
        csl = self%lib%slingo_adapted%tab(:,5)
        dsl = self%lib%slingo_adapted%tab(:,6)
        esl = self%lib%slingo_adapted%tab(:,3)
        fsl = self%lib%slingo_adapted%tab(:,4)

        ! Build daylight index
        nvalid = count(self%solz < 90.0d0)
        allocate(daylight_idx(nvalid))
        nvalid = 0
        do i = 1, self%p_size
            if (self%solz(i) < 90.0d0) then
                nvalid = nvalid + 1
                daylight_idx(nvalid) = i
            end if
        end do

#ifdef _OPENACC
            devtype = acc_device_nvidia
            ! ngpus = acc_get_num_devices(devtype)
            ngpus = 1
        if (ngpus < 1) ngpus = 1
#else
            ngpus = 1
#endif
        chunk = ceiling(real(nvalid) / ngpus)

        ! print *, "Using", ngpus, "GPUs with chunk size:", chunk

        ! Launch all GPU kernels asynchronously
        do g = 0, ngpus-1
#ifdef _OPENACC
            call acc_set_device_num(g, devtype)
#endif
            start_idx = g * chunk + 1
            end_idx = min((g + 1) * chunk, nvalid)
            local_chunk_size = end_idx - start_idx + 1
            
            if (start_idx > nvalid) cycle  ! Skip if no work for this GPU
            
            
            ! Allocate per-GPU working arrays
            allocate(ed_local(rows_size), es_local(rows_size))
            allocate(ta_local(rows_size), wa_local(rows_size), asym_local(rows_size))
            allocate(rlamu_local(rows_size), td_local(rows_size), ts_local(rows_size))
            allocate(tcd_local(rows_size), tcs_local(rows_size), tgas_local(rows_size))
            allocate(edclr_local(rows_size), esclr_local(rows_size))
            allocate(edcld_local(rows_size), escld_local(rows_size))
            
            ! Initialize from main arrays
            rlamu_local = self%rlamu
            td_local = self%td
            ts_local = self%ts
            tcd_local = self%tcd
            tcs_local = self%tcs
            tgas_local = self%tgas
            edclr_local = self%edclr
            esclr_local = self%esclr
            edcld_local = self%edcld
            escld_local = self%escld
#ifdef _OPENACC
            !$acc data create(ta_local, wa_local, asym_local, rlamu_local, &
            !$acc&            td_local, ts_local, tcd_local, tcs_local, tgas_local, &
            !$acc&            ed_local, es_local, edclr_local, esclr_local, &
            !$acc&            edcld_local, escld_local) &
            !$acc& copyin(daylight_idx(start_idx:end_idx), &
            !$acc&        oz, wv, rh, slp, wsm, ccov, rlwp, cdre, &
            !$acc&        asymp, ssalb, taua, &
            !$acc&        oza, ao, awv, dsl, csl, tab2, fobar, &
            !$acc&        bsl, esl, fsl, asl, aco2) &
            !$acc& async(g)
            !$acc parallel loop gang vector vector_length(64) async(g)
#else 
            !$OMP PARALLEL DO 
            !$OMP& SCHEDULE(DYNAMIC)
#endif
            do j = start_idx, end_idx
                i = daylight_idx(j)
                cosunz = cos(self%solz(i) * rad_1)
                sunz = self%solz(i)

                pres = slp(i)
                ws = wsm(i)
                ozone = oz(i)
                wvapor = wv(i)
                relhum = rh(i)
                
                ta_local(:) = taua(i,:)
                asym_local(:) = asymp(i,:)
                wa_local(:) = ssalb(i,:)
                rlamu_local(:) = self%rlamu(:)
                cov = ccov(i)
                clwp = rlwp(i)
                re = cdre(i)
                call light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
                        am, vi, cov, clwp, re, rows_size, &
                        fobar, oza, awv, ao, aco2, tab2, &
                        asl, bsl, csl, dsl, esl, fsl, &
                        ta_local, wa_local, asym_local, rlamu_local, &
                        td_local, ts_local, tcd_local, tcs_local, tgas_local, &
                        ed_local, es_local, edclr_local, esclr_local, &
                        edcld_local, escld_local, error)
                self%eda(i,:) = ed_local(:)
                self%esa(i,:) = es_local(:)
            end do
#ifdef _OPENACC
            !$acc end parallel loop
            !$acc update host(ed_local, es_local) async(g)
            !$acc end data
#else
            !$OMP END PARALLEL DO
#endif

            ! Deallocate per-GPU arrays after sync
            deallocate(ed_local, es_local)
            deallocate(ta_local, wa_local, asym_local)
            deallocate(rlamu_local, td_local, ts_local)
            deallocate(tcd_local, tcs_local, tgas_local)
            deallocate(edclr_local, esclr_local, edcld_local, escld_local)
        end do

        ! Wait for all GPUs to complete
#ifdef _OPENACC
        do g = 0, ngpus-1
            call acc_set_device_num(g, devtype)
            !$acc wait(g)
        end do
#endif
        
        ! Clean up shared arrays
        deallocate(fobar, oza, awv, ao, aco2, tab2)
        deallocate(asl, bsl, csl, dsl, esl, fsl)
        deallocate(daylight_idx)
    end subroutine sfcirr
end submodule oasim_sfcirr