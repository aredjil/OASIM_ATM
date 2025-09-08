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

        real(kind=real_kind), dimension(2) :: tdb_tdir
        real(kind=real_kind), dimension(self%lib%rows) :: input
        
        real(kind=real_kind), dimension(3, self%lib%rows):: input_2


        fobar => self%lib%atmo_adapted%tab(:,1)
        oza => self%lib%atmo_adapted%tab(:,3)
        awv => self%lib%atmo_adapted%tab(:,4)
        ao => self%lib%atmo_adapted%tab(:,5)
        aco2 => self%lib%atmo_adapted%tab(:,6)

        if (pres < 0.0e0 .or. ws < 0.0e0 .or. relhum < 0.0e0 &
        .or. ozone < 0.0e0 .or. wvapor < 0.0e0) then
            self%ed = 0
            self%es = 0
        end if

        rtmp = (93.885e0 - sunz) ** (-1.253e0)
        rmu0 = cosunz + 0.15e0 * rtmp
        rm = 1.0e0 / rmu0
        otmp = (cosunz * cosunz + ozfac1) ** 0.5e0
        rmo = ozfac2 / otmp

        rmp = pres / p0 * rm
        
        !$acc parallel loop 
        do concurrent (i = 1:self%lib%rows)
            to = oza(i) * ozone * 1.0e-3
            oarg = -to * rmo

            ag = ao(i) + aco2(i)
            gtmp = (1.0e0 + 118.3e0 * ag * rmp) ** 0.45e0
            gtmp2 = -1.41e0 * ag * rmp
            garg = gtmp2 / gtmp

            wtmp = (1.0e0 + 20.07e0 * awv(i) * wvapor * rm) ** 0.45e0
            wtmp2 = -0.2385e0 * awv(i) * wvapor * rm
            warg = wtmp2 / wtmp
            self%tgas(i) = exp(oarg + garg + warg)
        end do
        ! Here I think I should pass the attributes as arguments 
        ! to the ctrans function ? 
        input_2 = self%clrtrans(cosunz, rm, rmp, ws, relhum, am, vi)
        ! Then return the ouputs and assign them to their respective 
        ! Attributes td, ts and ta 
        self%td = input_2(1, :)
        self%ts = input_2(3, :)
        self%ta = input_2(3, :)

        self%edclr = daycor * cosunz * fobar * self%tgas * self%td
        self%esclr = daycor * cosunz * fobar * self%tgas * self%ts

        ! tdb_tdir here contain tdb and tdir 
        tdb_tdir = self%slingo(rmu0, clwp, re)


        self%tcd = tdb_tdir(1)
        self%tcs = tdb_tdir(2)        

        self%edcld = daycor * cosunz * fobar * self%tgas * self%tcd
        self%escld = daycor * cosunz * fobar * self%tgas * self%tcs

        ccov1 = cov * 1.0e-2
        ! The retun values of the light submodule ? 
        self%ed = (1.0e0 - ccov1) * self%edclr + ccov1 * self%edcld
        self%es = (1.0e0 - ccov1) * self%esclr + ccov1 * self%escld
    end subroutine light
end submodule