submodule (oasim) oasim_submod
    implicit none

contains
    module subroutine sunmod(self, iday, iyr, sec, points, rs)
    !!! NOTICE: radeg has been read from constants
    !!! NOTICE: here gmt has been replaced by sec
        implicit none

        class(calc_unit) :: self
        real(kind=real_kind), intent(in) :: sec
        integer, intent(in) :: iday, iyr
        integer, dimension(:), intent(in) :: points
        real(kind=real_kind), intent(out) :: rs

        real(kind=real_kind) :: day, gha, ghar
        real(kind=real_kind), dimension(3) :: suni, sung

        call sun2000(iyr, iday, sec, suni, rs)

        day = float(iday) + sec / 86400.0e0
        call gha2000(iyr, day, gha)
        ghar = gha * rad_1

        sung(1) = suni(1) * cos(ghar) + suni(2) * sin(ghar)
        sung(2) = suni(2) * cos(ghar) - suni(1) * sin(ghar)
        sung(3) = suni(3)

        self%sunv = matmul(self%lib%up(points,:), sung)
        self%sunn = matmul(self%lib%no(points,:), sung)
        self%sune = matmul(self%lib%ea(points,:), sung)

        self%sunz = rad * atan2(sqrt(self%sunn * self%sunn + self%sune * self%sune), self%sunv)
     end subroutine sunmod

    subroutine sun2000(iyr, iday, sec, sunvec, rs)
    !!! NOTICE: radeg has been read from constants
        implicit none

        real(kind=real_kind), intent(in) :: sec
        integer, intent(in) :: iyr, iday
        real(kind=real_kind), dimension(3), intent(out) :: sunvec
        real(kind=real_kind), intent(out) :: rs

        real(kind=real_kind), parameter :: xk = 5.6932e-3

        integer :: imon, nt, nutime
        real(kind=real_kind) :: rjd, t, xls, gs, xlm, omega, dpsi, eps 
        real(kind=real_kind) :: g2, g4, g5, dls, xlsg, xlsa

        imon = 1
        nutime = -99999

        rjd = float(jd(iyr, imon, iday))
        t = rjd - 2451545.0e0 + (sec - 43200.0e0) / 86400.0e0

        call ephparms(t, xls, gs, xlm, omega)
        nt = int(t)
        if (nt /= nutime) then
            nutime = nt
            call nutate(t, xls, gs, xlm, omega, dpsi, eps)
        end if

        g2 = 50.40828e0 + 1.60213022e0 * t
        g2 = mod(g2, 360.0e0)

        g4 = 19.38816e0 + 0.52402078e0 * t
        g4 = mod(g4, 360.0e0)
        
        g5 = 20.35116e0 + 0.08309121e0 * t
        g5 = mod(g5, 360.0e0)

        rs = +1.00014e0 - 0.01671e0 * cos(gs * rad_1) & 
             -0.00014e0 * cos(2.0e0 * gs * rad_1)

        dls = (6893.0e0 - 4.6543463e-4 * t) * sin(gs * rad_1) &
              +  72.0e0 * sin(2.0e0 * gs * rad_1) &
              -   7.0e0 * cos((gs - g5) * rad_1) &
              +   6.0e0 * sin((xlm - xls) * rad_1) &
              +   5.0e0 * sin((4.0e0 * gs - 8.0e0 * g4 + 3.0e0 * g5) * rad_1) &
              -   5.0e0 * cos((2.0e0 * gs - 2.0e0 * g2) * rad_1) &
              -   4.0e0 * sin((gs - g2) * rad_1) &
              +   4.0e0 * cos((4.0e0 * gs - 8.0e0 * g4 + 3.0e0 * g5) * rad_1) &
              +   3.0e0 * sin((2.0e0 * gs - 2.0e0 * g2) * rad_1) &
              -   3.0e0 * sin(g5 * rad_1) &
              -   3.0e0 * sin((2.0e0 * gs - 2.0e0 * g5) * rad_1)
              
        xlsg = xls + dls / 3600.0e0
        xlsa = xlsg + dpsi - xk / rs

        sunvec(1) = cos(xlsa * rad_1)
        sunvec(2) = sin(xlsa * rad_1) * cos(eps * rad_1)
        sunvec(3) = sin(xlsa * rad_1) * sin(eps * rad_1)
    end subroutine

    subroutine gha2000(iyr, day, gha)
    !!! NOTICE: radeg has been read from constants
        implicit none

        real(kind=real_kind), intent(in) :: day
        integer, intent(in) :: iyr
        real(kind=real_kind), intent(out) :: gha

        integer :: imon, nutime, iday, jday, nt
        real(kind=real_kind) :: fday, gmst, xls, gs, xlm, omega, dpsi, eps, t

        imon = 1
        nutime = -99999

        iday = int(day)
        fday = day - iday
        jday = jd(iyr, imon, iday)
        t = jday - 2451545.5e0 + fday

        gmst = 100.4606184e0 + 0.9856473663e0 * t + 2.908e-13 * t * t

        nt = int(t)
        if (nt /= nutime) then
            nutime = nt
            call ephparms(t, xls, gs, xlm, omega)
            call nutate(t, xls, gs, xlm, omega, dpsi, eps)
        end if

        gha = gmst + dpsi * cos(eps * rad_1) + fday * 360.0e0
        gha = mod(gha, 360.0e0)
        if (gha < 0.0e0) gha = gha + 360.0e0
    end subroutine gha2000

    pure integer function jd(y, m, d)
        implicit none

        integer, intent(in) :: y, m ,d

        jd = 367 * y - 7 * (y + (m + 9) / 12) / 4 + 275 * m / 9 + d + 1721014
    end function jd

    subroutine ephparms(t, xls, gs, xlm, omega)
        implicit none

        real(kind=real_kind), intent(in) :: t
        real(kind=real_kind), intent(out) :: xls, gs, xlm, omega

        xls = 280.46592e0 + 0.9856473516e0 * t
        xls = mod(xls, 360.0e0)

        gs = 357.52772e0 + 0.9856002831e0 * t
        gs = mod(gs, 360.0e0)
        
        xlm = 218.31643e0 + 13.17639648e0 * t 
        xlm = mod(xlm, 360.0e0)
        
        omega = 125.04452e0 - 0.0529537648e0 * t 
        omega = mod(omega, 360.0e0)        
    end subroutine ephparms

    subroutine nutate(t, xls, gs, xlm, omega, dpsi, eps)
        implicit none

        real(kind=real_kind), intent(in) :: t, xls, gs, xlm, omega
        real(kind=real_kind), intent(out) :: dpsi, eps

        real(kind=real_kind) :: epsm, deps

        dpsi = -17.1996e0 * sin(omega * rad_1) &
               +0.2062e0 * sin(2.0d0 * omega * rad_1) &
               -1.3187e0 * sin(2.0d0 * xls * rad_1) &
               +0.1426e0 * sin(gs * rad_1) &
               -0.2274e0 * sin(2.0d0 * xlm * rad_1)

        epsm = 23.439291e0 - 3.560e-7 * t
        deps = 9.2025e0 * cos(omega * rad_1) + 0.5736d0 * cos(2.0e0 * xls * rad_1)
        eps = epsm + deps / 3600.0e0
        dpsi = dpsi / 3600.0e0
    end subroutine nutate
end submodule oasim_submod