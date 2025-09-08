submodule (oasim) oasim_ocalbedo
    implicit none

contains
    module subroutine ocalbedo(self, wsm)
    !!! lib vars: rod, ros, wfac, solz
        implicit none

        class(calc_unit) :: self
        real(kind=real_kind), dimension(:), intent(in) :: wsm

        real(kind=real_kind), parameter :: rn = 1.341e0
        real(kind=real_kind), parameter :: roair = 1.2e3

        integer :: i
        real(kind=real_kind) :: cn, rof, rosps, rospd, rtheta, sintr, tanrmin, tanrpls, a, b
        real(kind=real_kind) :: rthetar, rmin, rpls, sinrmin, sinrpls, sinp, tanp

        do i = 1, self%p_size
            if (wsm(i) > 4.0e0) then
                if (wsm(i) <= 7.0e0) then
                    cn = 6.2e-4 + 1.56e-3 / wsm(i)
                    rof = roair * cn * 2.2d-5 * wsm(i) * wsm(i) - 4.0e-4
                else
                    cn = 0.49d-3 + 6.5d-5 * wsm(i)
                    rof = (roair * cn * 4.5e-5 - 4.0e-5) * wsm(i) * wsm(i)
                end if
                rosps = 5.7e-2
            else
                rof = 0.0d0
                rosps = 6.6e-2
            end if

            if (self%solz(i) < 40.0e0 .or. wsm(i) < 2.0e0) then
                if (self%solz(i) == 0.0e0) then  !!! put trhreshold !!!
                    rospd = 2.11e-2
                else
                    rtheta = self%solz(i) * rad_1
                    sintr = sin(rtheta) / rn
                    rthetar = asin(sintr)
                    rmin = rtheta - rthetar
                    rpls = rtheta + rthetar
                    sinrmin = sin(rmin)
                    sinrpls = sin(rpls)
                    tanrmin = tan(rmin)
                    tanrpls = tan(rpls)
                    sinp = (sinrmin * sinrmin) / (sinrpls * sinrpls)
                    tanp = (tanrmin * tanrmin) / (tanrpls * tanrpls)
                    rospd = 0.5d0 * (sinp + tanp)
                end if
            else
                a = 2.53e-2
                b = -7.14e-4 * wsm(i) + 6.18e-2
                rospd = a * exp(b * (self%solz(i) - 40.0e0))
            end if

            self%rod(i,:) = rospd + rof * self%wfac
            self%ros(i,:) = rosps + rof * self%wfac
        end do
    end subroutine ocalbedo
end submodule oasim_ocalbedo