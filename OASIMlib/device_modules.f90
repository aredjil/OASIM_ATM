module oasim_device
    use oasim_common, only: real_kind
    use, intrinsic:: iso_fortran_env, only: error_unit

    implicit none
contains
    !$acc routine (light) seq
subroutine light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
                            am, vi, cov, clwp, re, rows,                     &
                            fobar, oza, awv, ao, aco2,                       &
                            tab2, asl, bsl, csl, dsl, esl, fsl,              &
                            ta, wa, asym, rlamu,                             &
                            td, ts, tcd, tcs, tgas,                          &
                            ed, es, edclr, esclr, edcld, escld,              &
                            error)
    implicit none
    ! Inputs
    real(kind=real_kind), intent(in) :: sunz, cosunz, daycor, pres, ws, ozone
    real(kind=real_kind), intent(in) :: wvapor, relhum, am, vi, cov, clwp, re
    integer, intent(in)      :: rows
    real(kind=real_kind), intent(in) :: fobar(rows), oza(rows), awv(rows), ao(rows), aco2(rows)
    real(kind=real_kind), intent(in) :: tab2(rows)
    real(kind=real_kind), intent(in) :: asl(rows), bsl(rows), csl(rows), dsl(rows), esl(rows), fsl(rows)

    real(kind=real_kind), intent(inout) :: ta(rows), wa(rows), asym(rows)
    real(kind=real_kind), intent(in)    :: rlamu(rows)
    
    ! In/Out arrays
    real(kind=real_kind), intent(out) :: td(rows), ts(rows), tcd(rows), tcs(rows), tgas(rows)
    real(kind=real_kind), intent(out) :: ed(rows), es(rows), edclr(rows), esclr(rows), edcld(rows), escld(rows)

    ! Other outputs
    logical, intent(out) :: error

    ! Local variables
    real(kind=real_kind), parameter :: ozfac1 = 44.0d0 / 6370.d0
    real(kind=real_kind), parameter :: ozfac2 = 1.0d0 + 22.0d0 / 6370.0d0
    real(kind=real_kind), parameter :: p0 = 1013.25d0

    real(kind=real_kind) :: rtmp, rmu0, rm, otmp, rmo, rmp, to, oarg, ag
    real(kind=real_kind) :: gtmp, gtmp2, garg, wtmp, wtmp2, warg, ccov1
    integer :: i

    if (pres < 0.0d0 .or. ws < 0.0d0 .or. relhum < 0.0d0 .or. ozone < 0.0d0 .or. wvapor < 0.0d0) then
        ed = 0.0d0
        es = 0.0d0
        return
    end if

    rtmp = (93.885d0 - sunz) ** (-1.253d0)
    rmu0 = cosunz + 0.15d0 * rtmp
    rm = 1.0d0 / rmu0
    otmp = (cosunz * cosunz + ozfac1) ** 0.5d0
    rmo = ozfac2 / otmp
    rmp = pres / p0 * rm

    do i = 1, rows
        to = oza(i) * ozone * 1.0d-3
        oarg = -to * rmo

        ag = ao(i) + aco2(i)
        gtmp  = (1.0d0 + 118.3d0 * ag * rmp) ** 0.45d0
        gtmp2 = -1.41d0 * ag * rmp
        garg  = gtmp2 / gtmp

        wtmp  = (1.0d0 + 20.07d0 * awv(i) * wvapor * rm) ** 0.45d0
        wtmp2 = -0.2385d0 * awv(i) * wvapor * rm
        warg  = wtmp2 / wtmp

        tgas(i) = exp(oarg + garg + warg)
    end do

    ! Atmospheric transmission
    ! call clrtrans(cosunz, rm, rmp, ws, relhum, am, vi, rows, tab2, td, ts, error)
    call clrtrans(cosunz, rm, rmp, ws, relhum, am, vi, rows, tab2, ta, wa, asym, rlamu, td, ts, error)
    edclr = daycor * cosunz * fobar * tgas * td
    esclr = daycor * cosunz * fobar * tgas * ts

    ! Slingo parameterization
    call slingo(rmu0, clwp, re, rows, asl, bsl, csl, dsl, esl, fsl, tcd, tcs)

    edcld = daycor * cosunz * fobar * tgas * tcd
    escld = daycor * cosunz * fobar * tgas * tcs

    ccov1 = cov * 1.0d-2
    ed = (1.0d0 - ccov1) * edclr + ccov1 * edcld
    es = (1.0d0 - ccov1) * esclr + ccov1 * escld
end subroutine light

!$acc routine (slingo) seq
subroutine slingo(rmu0, clwp, cre, rows, asl, bsl, csl, dsl, esl, fsl, tcd, tcs)
    !!! Standalone version - all variables passed as parameters
        use :: oasim_common, only: real_kind
        implicit none

        ! Input parameters
        real(kind=real_kind), intent(in) :: rmu0, clwp, cre
        integer, intent(in) :: rows
        real(kind=real_kind), dimension(rows), intent(in) :: asl, bsl, csl, dsl, esl, fsl
        
        ! Output arrays
        real(kind=real_kind), dimension(rows), intent(out) :: tcd, tcs

        ! Local variables
        real(kind=real_kind) :: re, tauc, oneomega, omega, g, b0, bmu0, f, u2, sqarg
        real(kind=real_kind) :: eps, rm, e, val1, val2, val3, rnum, rden, gama1, gama2
        real(kind=real_kind) :: tdb, rdif, tdif, tdir
        real(kind=real_kind), dimension(4) :: alpha
        real(kind=real_kind), parameter :: const_3f7 = 3.0 / 7.0
        real(kind=real_kind), parameter :: const_7f4 = 7.0 / 4.0
        integer :: i

        re = (10.0 + 11.8) * 0.5
        if (cre >= 0.0) re = cre

        do i = 1, rows
            tauc = clwp * (asl(i) * 1.0d-2 + bsl(i) / re)
            oneomega = csl(i) + dsl(i) * re
            omega = 1.0d0 - oneomega
            g = esl(i) + fsl(i) * 1.0d-3 * re
            b0 = const_3f7 * (1.0d0 - g)
            bmu0 = 0.5 - 0.75 * rmu0 * g / (1.0d0 + g)
            f = g * g
            u2 = const_7f4 * (1.0d0 - ((1.0d0 - omega) / (7.0 * omega * b0)))
            u2 = max(u2, 0.0d0)
            alpha(1) = const_7f4 * (1.0d0 - omega * (1.0d0 - b0))
            alpha(2) = u2 * omega * b0
            alpha(3) = (1.0d0 - f) * omega * bmu0
            alpha(4) = (1.0d0 - f) * omega * (1.0d0 - bmu0)
            sqarg = alpha(1) * alpha(1) - alpha(2) * alpha(2)
            sqarg = max(sqarg, 1.0d-9)
            eps = sqrt(sqarg)
            rm = alpha(2) / (alpha(1) + eps)
            e = exp(-eps * tauc)
            val1 = 1.0d0 - omega * f
            val2 = eps * eps * rmu0 * rmu0
            rnum = val1 * alpha(3) - rmu0 * (alpha(1) * alpha(4) + alpha(2) * alpha(3))
            rden = val1 * val1 - val2
            gama1 = rnum / (rden + 1.0d-9)
            rnum = -val1 * alpha(4) - rmu0 * (alpha(1) * alpha(4) + alpha(2) * alpha(3))
            gama2 = rnum / (rden + 1.0d-9)
            tdb = exp(-val1 * tauc / rmu0)
            val3 = 1.0d0 - E * E * rm * rm
            rdif = rm * (1.0d0 - E * E) / (val3 + 1.0d-9)
            tdif = E * (1.0d0 - rm * rm) / (val3 + 1.0d-9)
            tdir = -gama2 * tdif - gama1 * tdb * rdif + gama2 * tdb
            tcd(i) = tdb
            tcs(i) = tdir
        end do

end subroutine slingo

    !$acc routine (clrtrans) seq
    subroutine clrtrans(cosunz, rm, rmp, ws, relhum, am, vi, &
                                         rows, thray, ta, wa, asym, rlamu, td, ts, error)
    !!! Standalone version - all variables passed as parameters
        implicit none

        ! Input parameters
        real(kind=real_kind), intent(in) :: cosunz, rm, rmp, ws, relhum, am, vi
        integer, intent(in) :: rows
        real(kind=real_kind), dimension(rows), intent(in) :: thray, rlamu
        
        ! Input/Output arrays
        real(kind=real_kind), dimension(rows), intent(inout) :: ta, wa, asym
        
        ! Output arrays
        real(kind=real_kind), dimension(rows), intent(out) :: td, ts
        logical, intent(out) :: error

        ! Local variables
        real(kind=real_kind) :: beta, eta, wa1, afs, bfs, rtra, omegaa, alg, fa
        real(kind=real_kind) :: tarm, atra, taa, tas, dray, daer
        integer :: i

        call navaer(relhum, am, vi, ws, beta, eta, wa1, afs, bfs)

        error = .false.
        do i = 1, rows
            rtra = exp(-thray(i) * rmp)
            if (ta(i) < 0.0d0) then
                ta(i) = beta * rlamu(i) ** eta
                ! perhaps one can bring it outside the loop and use where...
            end if

            if (wa(i) < 0.0d0) then
                omegaa = wa1
            else
                omegaa = wa(i)
            end if

            if (asym(i) >= 0.0d0) then
                alg = log(1.0d0 - asym(i) + 1.0d-2)
                afs = alg * (1.459d0 + alg * (0.1595d0 + alg * 0.4129d0))
                bfs = alg * (0.0783d0 + alg * (-0.3824d0 - alg * 0.5874d0))
            end if

            ! if (ta(i) < 0.0d0 .or. omegaa < 0.0d0) then
            !     write(error_unit, *) "ERROR in ta or omegaa"
            !     write(error_unit, *) "nl, ta, wa, asym = ", i, ta(i), wa(i), asym(i)
            !     error = .true.
            ! end if

            ! fa = 1.0d0 - 0.50d0 * exp((afs + bfs * cosunz) * cosunz)
            ! if (fa < 0.0d0) then
            !     write(error_unit, *) "ERROR in Fa"
            !     write(error_unit, *) "nl, ta, wa, asym = ", i, ta(i), wa(i), asym(i)
            !     error = .true.
            ! end if

            tarm = ta(i) * rm
            atra = exp(-tarm)
            taa = exp(-(1.0d0 - omegaa) * tarm)
            tas = exp(-omegaa * tarm)

            td(i) = rtra * atra

            dray = taa * 0.5d0 * (1.0d0 - rtra ** 0.95d0)
            daer = rtra ** 1.5d0 * taa * fa * (1.0d0 - tas)
            ts(i) = dray + daer
        end do
    end subroutine clrtrans

!$acc routine (navaer) seq
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
        
        do i = 1, 3
            dndr(i) = 0.0d0
            do j = 1, 3
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
        do i = 1, 3
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

end module oasim_device