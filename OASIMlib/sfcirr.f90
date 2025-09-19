submodule (oasim) oasim_sfcirr
    use oasim_device, only: light, clrtrans, slingo
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
        real(kind=real_kind), dimension(self%lib%rows) :: ed_local, es_local
        integer :: i, rows_size

        ! Local copies of frequently accessed derived type components
        real(kind=real_kind), dimension(:), allocatable :: fobar, oza, awv, ao, aco2, tab2
        real(kind=real_kind), dimension(:), allocatable :: asl, bsl, csl, dsl, esl, fsl
        real(kind=real_kind), dimension(:), allocatable :: ta_local, wa_local, asym_local
        real(kind=real_kind), dimension(:), allocatable :: rlamu_local, td_local, ts_local
        real(kind=real_kind), dimension(:), allocatable :: tcd_local, tcs_local, tgas_local

        am = self%lib%init_parameters%am
        vi = self%lib%init_parameters%vi
        rows_size = self%lib%rows
        
        ! Initialize output arrays
        self%eda = 0.0d0
        self%esa = 0.0d0
        
        rday = real(iday, real_kind) + sec_c * daypersec
        daycor = 1.0 + 1.67d-2 * cos(pi2 * (rday - 3.0d0) / 365.0d0)
        daycor = daycor * daycor

        ! Allocate and copy data to local arrays to avoid derived type issues
        allocate(fobar(rows_size), oza(rows_size), awv(rows_size), ao(rows_size))
        allocate(aco2(rows_size), tab2(rows_size))
        allocate(asl(rows_size), bsl(rows_size), csl(rows_size), dsl(rows_size))
        allocate(esl(rows_size), fsl(rows_size))
        allocate(ta_local(rows_size), wa_local(rows_size), asym_local(rows_size))
        allocate(rlamu_local(rows_size), td_local(rows_size), ts_local(rows_size))
        allocate(tcd_local(rows_size), tcs_local(rows_size), tgas_local(rows_size))

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

        ! Copy initial values for arrays that will be modified
        rlamu_local = self%rlamu
        td_local = self%td
        ts_local = self%ts
        tcd_local = self%tcd
        tcs_local = self%tcs
        tgas_local = self%tgas
        !$acc data create(ta_local, wa_local, asym_local, rlamu_local, td_local, ts_local, tcd_local, tcs_local, tgas_local, ed_local, es_local)
        !$acc parallel loop 
        do i = 1, self%p_size
            cosunz = cos(self%solz(i) * rad_1)
            sunz = self%solz(i)
            
            if (sunz < 90.0d0) then
                pres = slp(i)
                ws = wsm(i)
                ozone = oz(i)
                wvapor = wv(i)
                relhum = rh(i)
                
                ! Copy aerosol properties for this point
                ta_local(:) = taua(i,:)
                asym_local(:) = asymp(i,:)
                wa_local(:) = ssalb(i,:)
                
                cov = ccov(i)
                clwp = rlwp(i)
                re = cdre(i)

                call light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
                          am, vi, cov, clwp, re, rows_size, &
                          fobar, oza, awv, ao, aco2, tab2, &
                          asl, bsl, csl, dsl, esl, fsl, &
                          ta_local, wa_local, asym_local, rlamu_local, &
                          td_local, ts_local, tcd_local, tcs_local, tgas_local, &
                          ed_local, es_local, self%edclr, self%esclr, &
                          self%edcld, self%escld, error)
                
                self%eda(i,:) = ed_local(:)
                self%esa(i,:) = es_local(:)
            end if
        end do
        !$acc end parallel loop
        !$acc end data 
        ! Copy back modified arrays to derived type components
        self%rlamu = rlamu_local
        self%td = td_local
        self%ts = ts_local
        self%tcd = tcd_local
        self%tcs = tcs_local
        self%tgas = tgas_local

        ! Clean up
        deallocate(fobar, oza, awv, ao, aco2, tab2)
        deallocate(asl, bsl, csl, dsl, esl, fsl)
        deallocate(ta_local, wa_local, asym_local)
        deallocate(rlamu_local, td_local, ts_local, tcd_local, tcs_local, tgas_local)

    end subroutine sfcirr
end submodule oasim_sfcirr