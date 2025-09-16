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
        real(kind=real_kind) :: cov, clwp, re !, sirr
        real(kind=real_kind):: am, vi
        integer :: i !, j
        am = self%lib%init_parameters%am
        vi = self%lib%init_parameters%vi

        self%eda = 0.0d0
        self%esa = 0.0d0

        rday = real(iday, real_kind) + sec_c * daypersec
        daycor = 1.0 + 1.67d-2 * cos(pi2 * (rday - 3.0d0) / 365.0d0)
        daycor = daycor * daycor
        !$acc parallel loop gang vector present(self, slp, wsm, oz, wv, rh, taua, asymp, ssalb, ccov, rlwp, cdre)
        do i = 1, self%p_size
            cosunz = cos(self%solz(i) * rad_1)
            sunz = self%solz(i)

            if (sunz < 90.0d0) then
                pres = slp(i)
                ws = wsm(i)
                ozone = oz(i)
                wvapor = wv(i)
                relhum = rh(i)

                self%ta = taua(i,:)
                self%asym = asymp(i,:)
                self%wa = ssalb(i,:)

                cov = ccov(i)
                clwp = rlwp(i)
                re = cdre(i)

                ! call self%light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
                !                 am, vi, &
                !                 cov, clwp, re, error)
                call light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
                                      am, vi, cov, clwp, re, self%lib%rows,                 &
                                      self%lib%atmo_adapted%tab(:,1), &   ! fobar
                                      self%lib%atmo_adapted%tab(:,3), &   ! oza
                                      self%lib%atmo_adapted%tab(:,4), &   ! awv
                                      self%lib%atmo_adapted%tab(:,5), &   ! ao
                                      self%lib%atmo_adapted%tab(:,6), &   ! aco2
                                      self%lib%atmo_adapted%tab(:,2), &   ! tab2
                                      self%lib%slingo_adapted%tab(:,1), & ! asl
                                      self%lib%slingo_adapted%tab(:,2), & ! bsl
                                      self%lib%slingo_adapted%tab(:,5), & ! csl
                                      self%lib%slingo_adapted%tab(:,6), & ! dsl
                                      self%lib%slingo_adapted%tab(:,3), & ! esl
                                      self%lib%slingo_adapted%tab(:,4), & ! fsl
                                       self%ta, self%wa, self%asym, self%rlamu,            &
                                      self%td, self%ts, self%tcd, self%tcs, self%tgas,     &
                                      self%ed, self%es, self%edclr, self%esclr, self%edcld, self%escld, &
                                      error)
                ! sirr = 0.0
                self%eda(i,:) = self%ed
                self%esa(i,:) = self%es

                ! sirr = sum(self%eda(i,:)) + sum(self%esa(i,:))
            else
                ! sirr = 0.0
                self%eda(i,:) = 0.0
                self%esa(i,:) = 0.0
            end if
        end do
        !$acc end parallel loop
    end subroutine sfcirr
end submodule oasim_sfcirr