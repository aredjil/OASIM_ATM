submodule (oasim) oasim_light
    implicit none

contains

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

end submodule
