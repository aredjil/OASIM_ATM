submodule (oasim) oasim_clrtrans
    use, intrinsic:: iso_fortran_env, only: error_unit
    implicit none

contains
    module subroutine clrtrans(self, cosunz, rm, rmp, ws, relhum, am, vi, error)
    !!! lib vars: ta, wa, asym, td, ts, rlamu
        implicit none

        class(calc_unit) :: self
        real(kind=real_kind), intent(in) :: cosunz, rm, rmp, ws, relhum, am, vi
        logical, intent(out) :: error

        real(kind=real_kind) :: beta, eta, wa1, afs, bfs, rtra, omegaa, alg, fa
        real(kind=real_kind) :: tarm, atra, taa, tas, dray, daer
        real(kind=real_kind), dimension(:), pointer :: thray
        integer :: i

        thray => self%lib%atmo_adapted%tab(:,2)

        call navaer(relhum, am, vi, ws, beta, eta, wa1, afs, bfs)

        error = .false.

        do concurrent (i = 1:self%lib%rows)
            rtra = exp(-thray(i) * rmp)
            if (self%ta(i) < 0.0d0) then
                self%ta(i) = beta * self%rlamu(i) ** eta
                ! perhaps one can bring it outside the loop and use where...
            end if

            if (self%wa(i) < 0.0d0) then
                omegaa = wa1
            else
                omegaa = self%wa(i)
            end if

            if (self%asym(i) >= 0.0d0) then
                alg = log(1.0d0 - self%asym(i) + 1.0d-2)
                afs = alg * (1.459d0 + alg * (0.1595d0 + alg * 0.4129d0))
                bfs = alg * (0.0783d0 + alg * (-0.3824d0 - alg * 0.5874d0))
            end if

            if (self%ta(i) < 0.0d0 .or. omegaa < 0.0d0) then
                write(error_unit, *) "ERROR in ta or omegaa"
                write(error_unit, *) "nl, ta, wa, asym = ", i, self%ta(i), self%wa(i), self%asym(i)
                error = .true.
            end if

            fa = 1.0d0 - 0.50d0 * exp((afs + bfs * cosunz) * cosunz)
            if (fa < 0.0d0) then
                write(error_unit, *) "ERROR in Fa"
                write(error_unit, *) "nl, ta, wa, asym = ", i, self%ta(i), self%wa(i), self%asym(i)
                error = .true.
            end if

            tarm = self%ta(i) * rm
            atra = exp(-tarm)
            taa = exp(-(1.0d0 - omegaa) * tarm)
            tas = exp(-omegaa * tarm)

            self%td(i) = rtra * atra

            dray = taa * 0.5d0 * (1.0d0 - rtra ** 0.95d0)
            daer = rtra ** 1.5d0 * taa * fa * (1.0d0 - tas)
            self%ts(i) = dray + daer
        end do

    end subroutine clrtrans

    subroutine navaer(relhum, am, vi, ws, beta, eta, wa, afs, bfs)
        implicit none

        real(kind=real_kind), intent(in) :: relhum
        real(kind=real_kind), intent(in) :: am, vi, ws
        real(kind=real_kind), intent(out) :: beta, eta, wa, afs, bfs

        real(kind=real_kind), parameter :: rlam = 0.55d0
        real(kind=real_kind), dimension(3), parameter :: ro = [0.03d0, 0.24d0, 2.0d0]
        real(kind=real_kind), dimension(3), parameter :: r = [0.1d0, 1.0d0, 10.0d0]

        real(kind=real_kind), dimension(3) :: a, dndr
        real(kind=real_kind) :: rnum, rden, frh, arg, sumx, sumy, sumxy, sumx2, relhumnorm
        real(kind=real_kind) :: rlrn, rldndr, gama, alpha, rlogc, cext, asymp, alg, rval
        integer :: i, j

        relhumnorm = min(99.9d0, relhum)
        rnum = 2.0d0 - relhumnorm / 100.0d0
        rden = 6.0d0 * (1.0d0 - relhumnorm / 100.0d0)
        frh = (rnum / rden) ** 0.333

        a(1) = 2000.0d0 * am * am
        a(2) = 5.866d0 * (ws - 2.2d0)
        a(2) = max(0.5d0, a(2))
        a(3) = 0.01527d0 * (ws - 2.2d0) * 0.05d0
        a(3) = max(1.4d-5, a(3))

        do concurrent (i = 1:3)
            dndr(i) = 0.0d0
            do concurrent (j = 1:3)
                rden = frh * ro(j)
                arg = log(r(i) / rden)
                arg = arg * arg
                rval = a(j) * exp(-arg) / frh
                dndr(i) = dndr(i) + rval
            end do
        end do

        sumx = 0.0d0
        sumy = 0.0d0
        sumxy = 0.0d0
        sumx2 = 0.0d0
        do concurrent(i = 1:3)
            rlrn = log10(r(i))
            rldndr = log10(dndr(i))
            sumx = sumx + rlrn
            sumy = sumy + rldndr
            sumxy = sumxy + rlrn * rldndr
            sumx2 = sumx2 + rlrn * rlrn
        end do
        gama = sumxy / sumx2
        rlogc = sumy / 3.0d0 - gama * sumx / 3.0d0
        alpha = -(gama + 3.0d0)
        eta = -alpha

        cext = 3.91d0 / vi
        beta = cext * rlam ** alpha

        if (alpha > 1.2d0) then
            asymp = 0.65d0
        else if (alpha < 0.0d0) then
            asymp = 0.82d0
        else
            asymp = -0.14167d0 * alpha + 0.82d0
        end if

        alg = log(1.0d0 - asymp)
        afs = alg * (1.459d0 + alg * (0.1595d0 + alg * 0.4129d0))
        bfs = alg * (0.0783d0 + alg * (-0.3824d0 - alg * 0.5874d0))

        wa = (-3.2d-3 * am + 0.972d0) * exp(3.06d-4 * relhumnorm)
    end subroutine navaer
end submodule oasim_clrtrans
