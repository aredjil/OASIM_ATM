submodule (oasim) oasim_monrad
    implicit none

contains
    module subroutine monrad(self, points, iyr, iday, sec_b, sec_e, sp, msl, ws10, tco3, t2m, d2m, &
                             tcc, tclw, cdrem, taua, asymp, ssalb, edout, esout, error)

        implicit none

        class(calc_unit) :: self
        integer, dimension(:), intent(in) :: points
        integer, intent(in) :: iyr, iday
        real(kind=real_kind), intent(in) :: sec_b, sec_e
        real(kind=real_kind), dimension(:), intent(in) :: sp, msl, ws10, tco3, t2m, d2m, tcc, tclw, cdrem
        real(kind=real_kind), dimension(:, :), intent(in) :: taua, asymp, ssalb
        real(kind=real_kind), dimension(: ,:), intent(out) :: edout, esout
        logical, intent(out) :: error

        real(kind=real_kind), parameter :: a3 = 17.502d0 ! dimensionless
        real(kind=real_kind), parameter :: a4 = 32.19d0  ! Kelvin
        real(kind=real_kind), parameter :: To = 273.16d0 ! Kelvin
        real(kind=real_kind), parameter :: b1 = 0.14d0 * 0.01d0 * 611.21d0 ! cm
        real(kind=real_kind), parameter :: b2 = 0.21d0   ! cm

        real(kind=real_kind) :: sec_c

        sec_c = (sec_b + sec_e) * 0.5d0

        error = .true.

        if (size(points) /= self%p_size) call argument_error("points", size(points), "p_size", self%p_size)    
        if (size(sp) /= self%p_size) call argument_error("sp", size(sp), "p_size", self%p_size)
        if (size(ws10) /= self%p_size) call argument_error("ws10", size(ws10), "p_size", self%p_size)
        if (size(tco3) /= self%p_size) call argument_error("tco3", size(tco3), "p_size", self%p_size)
        if (size(t2m) /= self%p_size) call argument_error("t2m", size(t2m), "p_size", self%p_size)
        if (size(d2m) /= self%p_size) call argument_error("d2m", size(d2m), "p_size", self%p_size)
        if (size(tcc) /= self%p_size) call argument_error("tcc", size(tcc), "p_size", self%p_size)
        if (size(tclw) /= self%p_size) call argument_error("tclw", size(tclw), "p_size", self%p_size)
        if (size(cdrem) /= self%p_size) call argument_error("cdrem", size(cdrem), "p_size", self%p_size)
        if (size(taua, 1) /= self%p_size) call argument_error("taua, 1", size(taua, 1), "p_size", self%p_size)
        if (size(asymp, 1) /= self%p_size) call argument_error("asymp, 1", size(asymp, 1), "p_size", self%p_size)
        if (size(ssalb, 1) /= self%p_size) call argument_error("ssalb, 1", size(ssalb, 1), "p_size", self%p_size)
        if (size(edout, 1) /= self%p_size) call argument_error("edout, 1", size(edout, 1), "p_size", self%p_size)
        if (size(esout, 1) /= self%p_size) call argument_error("esout, 1", size(esout, 1), "p_size", self%p_size)
        if (size(taua, 2) /= self%lib%rows) call argument_error("taua, 2", size(taua, 2), "rows", self%lib%rows)
        if (size(asymp, 2) /= self%lib%rows) call argument_error("asymp, 2", size(asymp, 2), "rows", self%lib%rows)
        if (size(ssalb, 2) /= self%lib%rows) call argument_error("ssalb, 2", size(ssalb, 2), "rows", self%lib%rows)
        if (size(edout, 2) /= self%lib%rows) call argument_error("edout, 2", size(edout, 2), "rows", self%lib%rows)
        if (size(esout, 2) /= self%lib%rows) call argument_error("esout, 2", size(esout, 2), "rows", self%lib%rows)

        error = .false.  
        
        self%ed_b = exp(a3 * (d2m - To)/(d2m - a4))
        self%et_b = exp(-a3 * (t2m - To)/(t2m - a4))
        self%rh_b = 1.0d2 * self%ed_b * self%et_b
        self%wv_b = b1 * self%ed_b * sp / msl + b2


        call self%sfcsolz(iyr, iday, sec_b, sec_e, points)
        call self%ocalbedo(ws10)
        !!! after debug absorb numerical constant in the subroutine !!!
        call self%sfcirr(iday, sec_c, 1.0d-2 * sp, ws10, tco3, self%wv_b, self%rh_b, &
                         taua, asymp, ssalb, tcc, 1.0d3 * tclw, cdrem, error) 

        edout = self%eda * (1.0d0 - self%rod)
        esout = self%esa * (1.0d0 - self%ros)
    end subroutine monrad

    module subroutine monrad_debug(self, points, iyr, iday, sec_b, sec_e, slp, wsm, oz, wv, &
                                   rh, ccov, rlwp, cdre, taua, asymp, ssalb, edout, esout, error)
        implicit none

        class(calc_unit) :: self
        integer, dimension(:), intent(in) :: points
        integer, intent(in) :: iyr, iday
        real(kind=real_kind), intent(in) :: sec_b, sec_e
        real(kind=real_kind), dimension(:), intent(in) :: slp, wsm, oz, wv, rh, ccov, rlwp, cdre
        real(kind=real_kind), dimension(:, :), intent(in) :: taua, asymp, ssalb
        real(kind=real_kind), dimension(:, :), intent(out) :: edout, esout
        logical, intent(out) :: error

        real(kind=real_kind) :: sec_c
        
        !!! DEBUG SECTION !!!

        ! print *, "monrad_debug, input data"
        ! print *, "points = ", points
        ! print *, "year = ", iyr
        ! print *, "day = ", iday
        ! print *, "sec_b = ", sec_b
        ! print *, "sec_e = ", sec_e
        ! print *, "slp = ", slp
        ! print *, "wsp = ", wsm
        ! print *, "oz = ", oz
        ! print *, "wv = ", wv
        ! print *, "rh = ", rh
        ! print *, "ccov = ", ccov
        ! print *, "rlwp = ", rlwp
        ! print *, "cfre = ", cdre
        ! print "(*(G0,:,', '))", "taua = ", taua
        ! print "(*(G0,:,', '))", "asymp = ", asymp
        ! print "(*(G0,:,', '))", "ssalb = ", ssalb
        ! print *, "monrad_debug, end input data"
        
        
        sec_c = (sec_b + sec_e) * 0.5

        !!! Check Section (it may be commented out) !!!

        error = .true.

        if (size(points) /= self%p_size) return
        if (size(slp) /= self%p_size) return
        if (size(wsm) /= self%p_size) return
        if (size(oz) /= self%p_size) return
        if (size(wv) /= self%p_size) return
        if (size(rh) /= self%p_size) return
        if (size(ccov) /= self%p_size) return
        if (size(rlwp) /= self%p_size) return
        if (size(cdre) /= self%p_size) return
        if (size(taua, 1) /= self%p_size) return
        if (size(asymp, 1) /= self%p_size) return
        if (size(ssalb, 1) /= self%p_size) return
        if (size(edout, 1) /= self%p_size) return
        if (size(esout, 1) /= self%p_size) return
        if (size(taua, 2) /= self%lib%rows) return
        if (size(asymp, 2) /= self%lib%rows) return
        if (size(ssalb, 2) /= self%lib%rows) return
        if (size(edout, 2) /= self%lib%rows) return
        if (size(esout, 2) /= self%lib%rows) return

        error = .false.

        !!! Calculation Section !!!

        call self%sfcsolz(iyr, iday, sec_b, sec_e, points)
        call self%ocalbedo(wsm)
        !$acc enter data copyin(self)
        call self%sfcirr(iday, sec_c, slp, wsm, oz, wv, rh, taua, asymp, ssalb, ccov, rlwp, cdre, error)
        !$acc exit data delete(self)
        edout = self%eda * (1.0d0 - self%rod)
        esout = self%esa * (1.0d0 - self%ros)
    end subroutine monrad_debug

    subroutine argument_error(msg1, size1, msg2, size2)
        use, intrinsic:: iso_fortran_env, only: error_unit
        implicit none

        character(len=*), intent(in) :: msg1, msg2
        integer, intent(in) :: size1, size2

        write(error_unit, *) "ERROR: "//msg1//"=", size1, msg2//"=", size2
        stop
    end subroutine argument_error
end submodule oasim_monrad

