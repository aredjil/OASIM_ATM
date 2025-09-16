submodule (oasim) oasim_light
    implicit none

contains
    module subroutine light(self, sunz, cosunz, daycor, pres, ws, ozone, wvapor, &
        relhum, am, vi, cov, clwp, re, error)
        !!! lib vars: td, ts, tcd, tcs, ed, es, tgas
        !!!           edclr, esclr, edcld, escld, oza
        implicit none

        class(calc_unit) :: self
        real(kind=real_kind), intent(in) :: sunz, cosunz, daycor, pres, ws, ozone
        real(kind=real_kind), intent(in) :: wvapor, relhum, am, vi, cov, clwp, re
        logical, intent(out) :: error

        real(kind=real_kind), parameter :: ozfac1 = 44.0d0 / 6370.d0
        real(kind=real_kind), parameter :: ozfac2 = 1.0d0 + 22.0d0 / 6370.0d0
        real(kind=real_kind), parameter :: p0 = 1013.25d0

        real(kind=real_kind) :: rtmp, rmu0, rm, otmp, rmo, rmp, to, oarg, ag
        real(kind=real_kind) :: gtmp, gtmp2, garg, wtmp, wtmp2, warg, ccov1
        real(kind=real_kind), dimension(:), pointer :: fobar, oza, awv, ao, aco2
        integer :: i

        fobar => self%lib%atmo_adapted%tab(:,1)
        oza => self%lib%atmo_adapted%tab(:,3)
        awv => self%lib%atmo_adapted%tab(:,4)
        ao => self%lib%atmo_adapted%tab(:,5)
        aco2 => self%lib%atmo_adapted%tab(:,6)

        if (pres < 0.0d0 .or. ws < 0.0d0 .or. relhum < 0.0d0 &
        .or. ozone < 0.0d0 .or. wvapor < 0.0d0) then
            self%ed = 0
            self%es = 0
        end if

        rtmp = (93.885d0 - sunz) ** (-1.253d0)
        rmu0 = cosunz + 0.15d0 * rtmp
        rm = 1.0d0 / rmu0
        otmp = (cosunz * cosunz + ozfac1) ** 0.5d0
        rmo = ozfac2 / otmp

        rmp = pres / p0 * rm

        do i = 1, self%lib%rows
            to = oza(i) * ozone * 1.0d-3
            oarg = -to * rmo

            ag = ao(i) + aco2(i)
            gtmp = (1.0d0 + 118.3d0 * ag * rmp) ** 0.45d0
            gtmp2 = -1.41d0 * ag * rmp
            garg = gtmp2 / gtmp

            wtmp = (1.0d0 + 20.07d0 * awv(i) * wvapor * rm) ** 0.45d0
            wtmp2 = -0.2385d0 * awv(i) * wvapor * rm
            warg = wtmp2 / wtmp
            self%tgas(i) = exp(oarg + garg + warg)
        end do
        ! call self%clrtrans(cosunz, rm, rmp, ws, relhum, am, vi, error)
        call clrtrans(cosunz, rm, rmp, ws, relhum, am, vi, &
             self%lib%rows, &
             self%lib%atmo_adapted%tab(:,2), &
             self%ta, self%wa, self%asym, self%rlamu, &
             self%td, self%ts, error)
        self%edclr = daycor * cosunz * fobar * self%tgas * self%td
        self%esclr = daycor * cosunz * fobar * self%tgas * self%ts

        ! call self%slingo(rmu0, clwp, re)
        call slingo(rmu0, clwp, re, &
           self%lib%rows, &
           self%lib%slingo_adapted%tab(:,1), &  ! asl
           self%lib%slingo_adapted%tab(:,2), &  ! bsl
           self%lib%slingo_adapted%tab(:,5), &  ! csl
           self%lib%slingo_adapted%tab(:,6), &  ! dsl
           self%lib%slingo_adapted%tab(:,3), &  ! esl
           self%lib%slingo_adapted%tab(:,4), &  ! fsl
           self%tcd, self%tcs)

        self%edcld = daycor * cosunz * fobar * self%tgas * self%tcd
        self%escld = daycor * cosunz * fobar * self%tgas * self%tcs

        ccov1 = cov * 1.0d-2
        self%ed = (1.0d0 - ccov1) * self%edclr + ccov1 * self%edcld
        self%es = (1.0d0 - ccov1) * self%esclr + ccov1 * self%escld
    end subroutine light
end submodule