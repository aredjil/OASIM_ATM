submodule (oasim) oasim_slingo
    implicit none

contains
    module subroutine slingo(self, rmu0, clwp, cre)  
    !!! lib vars: tcd, tcs
        implicit none
        
        class(calc_unit) :: self
        real(kind=real_kind), intent(in) :: rmu0, clwp, cre

        real(kind=real_kind) :: re, tauc, oneomega, omega, g, b0, bmu0, f, u2, sqarg
        real(kind=real_kind) :: eps, rm, e, val1, val2, val3, rnum, rden, gama1, gama2
        real(kind=real_kind) :: tdb, rdif, tdif, tdir
        real(kind=real_kind), dimension(4) :: alpha
        real(kind=real_kind), dimension(:), pointer :: asl, bsl, csl, dsl, esl, fsl
        real(kind=real_kind), parameter :: const_3f7 = 3.0 / 7.0
        real(kind=real_kind), parameter :: const_7f4 = 7.0 / 4.0
        integer :: i

        asl => self%lib%slingo_adapted%tab(:,1)
        bsl => self%lib%slingo_adapted%tab(:,2)
        csl => self%lib%slingo_adapted%tab(:,5)
        dsl => self%lib%slingo_adapted%tab(:,6)
        esl => self%lib%slingo_adapted%tab(:,3)
        fsl => self%lib%slingo_adapted%tab(:,4)

        re = (10.0 + 11.8) * 0.5
        if (cre >= 0.0) re = cre
        !$acc parallel loop 
        do concurrent (i = 1:self%lib%rows)
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
            self%tcd(i) = tdb
            self%tcs(i) = tdir
        end do
    end subroutine slingo
end submodule oasim_slingo
