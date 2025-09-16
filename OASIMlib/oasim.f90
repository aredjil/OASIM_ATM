module oasim
        use :: configuration
        use :: tables
        use :: oasim_common, only: real_kind, string_length

        implicit none

        private

        type oasim_lib
            type(parameters) init_parameters
            type(table) atmo_adapted, abso_adapted, slingo_adapted
            integer :: rows
            real(kind=real_kind), dimension(:), pointer :: lat, lon
            integer :: m_size ! size of the mesh
            real(kind=real_kind), dimension(:,:), allocatable :: up, no, ea
        contains
            private

            procedure :: localvec
            procedure, public :: finalize => finalize_lib 
        end type oasim_lib

        interface oasim_lib
            module procedure :: init_lib
        end interface oasim_lib

        type calc_unit
            private

            type(oasim_lib), pointer :: lib
            integer :: p_size
            real(kind=real_kind), dimension(:), allocatable :: tcd, tcs, ta, asym, wa, rlamu, td, ts, wfac
            real(kind=real_kind), dimension(:), allocatable :: ed, es, tgas, edclr, esclr, edcld, escld
            real(kind=real_kind), dimension(:), allocatable :: solz, sunz, integ, sunv, sunn, sune
            real(kind=real_kind), dimension(:), allocatable :: ed_b, et_b, wv_b, rh_b
            real(kind=real_kind), dimension(:, :), allocatable :: rod, ros, eda, esa
        contains
            private

            procedure, public :: monrad
            procedure, public :: monrad_debug
            ! procedure :: slingo
            ! procedure :: clrtrans
            ! procedure :: light
            procedure :: ocalbedo
            procedure :: sunmod
            procedure :: sfcirr
            procedure :: sfcsolz
            procedure, public :: finalize => finalize_calc
        end type calc_unit

        interface calc_unit
            module procedure :: init_calc
        end interface calc_unit

        real(kind=real_kind), parameter :: pi = acos(-1.0d0)
        real(kind=real_kind), parameter :: pi2 = 2 * pi
        real(kind=real_kind), parameter :: rad = 180.0d0 / pi
        real(kind=real_kind), parameter :: rad_1 = pi / 180.0d0 

        ! interface
        !     module subroutine light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
        !                         am, vi, cov, clwp, re, rows,                     &
        !                         fobar, oza, awv, ao, aco2,                       &
        !                         tab2, asl, bsl, csl, dsl, esl, fsl,              &
        !                         ta, wa, asym, rlamu,                             &
        !                         td, ts, tcd, tcs, tgas,                          &
        !                         ed, es, edclr, esclr, edcld, escld,              &
        !                         error)
        !                 implicit none
        !                 real(kind=real_kind), intent(in) :: sunz, cosunz, daycor, pres, ws, ozone
        !                 real(kind=real_kind), intent(in) :: wvapor, relhum, am, vi, cov, clwp, re
        !                 integer, intent(in)      :: rows
        !                 real(kind=real_kind), intent(in) :: fobar(rows), oza(rows), awv(rows), ao(rows), aco2(rows)
        !                 real(kind=real_kind), intent(in) :: tab2(rows)
        !                 real(kind=real_kind), intent(in) :: asl(rows), bsl(rows), csl(rows), dsl(rows), esl(rows), fsl(rows)
                        
        !                 real(kind=real_kind), intent(inout) :: ta(rows), wa(rows), asym(rows)
        !                 real(kind=real_kind), intent(in)    :: rlamu(rows)
                        
        !                 real(kind=real_kind), intent(out) :: td(rows), ts(rows), tcd(rows), tcs(rows), tgas(rows)
        !                 real(kind=real_kind), intent(out) :: ed(rows), es(rows), edclr(rows), esclr(rows), edcld(rows), escld(rows)

        !                 logical, intent(out) :: error
        !      end subroutine light
        ! end interface

        interface
            ! module subroutine slingo(self, rmu0, clwp, cre)   
            !     use :: oasim_common, only: real_kind    
            !     class(calc_unit) :: self
            !     real(kind=real_kind), intent(in) :: rmu0, clwp, cre
            ! end subroutine slingo
            
            ! module subroutine slingo(rmu0, clwp, cre, rows, asl, bsl, csl, dsl, esl, fsl, tcd, tcs)
            !     use :: oasim_common, only: real_kind
            !     real(kind=real_kind), intent(in) :: rmu0, clwp, cre
            !     integer, intent(in) :: rows
            !     real(kind=real_kind), dimension(rows), intent(in) :: asl, bsl, csl, dsl, esl, fsl
            !     real(kind=real_kind), dimension(rows), intent(out) :: tcd, tcs
            ! end subroutine slingo

            ! module subroutine clrtrans(self, cosunz, rm, rmp, ws, relhum, am, vi, error)
            !     use :: oasim_common, only: real_kind 
            !     class(calc_unit) :: self
            !     real(kind=real_kind), intent(in) :: cosunz, rm, rmp, ws, relhum, am, vi
            !     logical, intent(out) :: error
            ! end subroutine clrtrans

            ! module subroutine clrtrans(cosunz, rm, rmp, ws, relhum, am, vi, &
            !               rows, thray, ta, wa, asym, rlamu, td, ts, error)
            !     use :: oasim_common, only: real_kind
            !     real(kind=real_kind), intent(in) :: cosunz, rm, rmp, ws, relhum, am, vi
            !     integer, intent(in) :: rows
            !     real(kind=real_kind), dimension(rows), intent(in) :: thray, rlamu
            !     real(kind=real_kind), dimension(rows), intent(inout) :: ta, wa, asym
            !     real(kind=real_kind), dimension(rows), intent(out) :: td, ts
            !     logical, intent(out) :: error
            ! end subroutine clrtrans


            ! module subroutine light(self, sunz, cosunz, daycor, pres, ws, ozone, wvapor, &
            !     relhum, am, vi, cov, clwp, re, error)
            !     use :: oasim_common, only: real_kind
            !     class(calc_unit) :: self
            !     real(kind=real_kind), intent(in) :: sunz, cosunz, daycor, pres, ws, ozone
            !     real(kind=real_kind), intent(in) :: wvapor, relhum, am, vi, cov, clwp, re              
            !     logical, intent(out) :: error
            ! end subroutine light

            module subroutine ocalbedo(self, wsm)
                use :: oasim_common, only: real_kind
                class(calc_unit) :: self
                real(kind=real_kind), dimension(:), intent(in) :: wsm
            end subroutine ocalbedo

            module subroutine sunmod(self, iday, iyr, sec, points, rs)
                use :: oasim_common, only: real_kind
                class(calc_unit) :: self
                real(kind=real_kind), intent(in) :: sec
                integer, intent(in) :: iday, iyr
                integer, dimension(:), intent(in) :: points
                real(kind=real_kind), intent(out) :: rs
            end subroutine sunmod
            !$acc routine seq
            module subroutine sfcirr(self, iday, sec_c, slp, wsm, oz, wv, rh, &
                taua, asymp, ssalb, ccov, rlwp, cdre, error)
                use :: oasim_common, only: real_kind     
                class(calc_unit) :: self
                integer, intent(in) :: iday                
                real(kind=real_kind), intent(in) :: sec_c
                real(kind=real_kind), dimension(:), intent(in) :: slp, wsm, oz, wv, rh
                real(kind=real_kind), dimension(:,:), intent(in) :: taua, asymp, ssalb
                real(kind=real_kind), dimension(:), intent(in) :: ccov, rlwp, cdre
                logical, intent(out) :: error
            end subroutine

            module subroutine sfcsolz(self, iyr, iday, sec_b, sec_e, points)
                use :: oasim_common, only: real_kind    
                class(calc_unit) :: self
                integer, intent(in) :: iyr, iday
                real(kind=real_kind), intent(in) :: sec_b, sec_e
                integer, dimension(:), intent(in) :: points
            end subroutine

            module subroutine monrad_debug(self, points, iyr, iday, sec_b, sec_e, slp, wsm, oz, wv, &
                rh, ccov, rlwp, cdre, taua, asymp, ssalb, edout, esout, error)
                use :: oasim_common, only: real_kind
                class(calc_unit) :: self
                integer, dimension(:), intent(in) :: points
                integer, intent(in) :: iyr, iday
                real(kind=real_kind), intent(in) :: sec_b, sec_e
                real(kind=real_kind), dimension(:), intent(in) :: slp, wsm, oz, wv, rh, ccov, rlwp, cdre
                real(kind=real_kind), dimension(:, :), intent(in) :: taua, asymp, ssalb
                real(kind=real_kind), dimension(:, :), intent(out) :: edout, esout
                logical, intent(out) :: error
            end subroutine monrad_debug

            module subroutine monrad(self, points, iyr, iday, sec_b, sec_e, sp, msl, ws10, tco3, t2m, d2m, &
                tcc, tclw, cdrem, taua, asymp, ssalb, edout, esout, error)
                use :: oasim_common, only: real_kind
                class(calc_unit) :: self
                integer, dimension(:), intent(in) :: points
                integer, intent(in) :: iyr, iday
                real(kind=real_kind), intent(in) :: sec_b, sec_e
                real(kind=real_kind), dimension(:), intent(in) :: sp, msl, ws10, tco3, t2m, d2m, tcc, tclw, cdrem
                real(kind=real_kind), dimension(:, :), intent(in) :: taua, asymp, ssalb
                real(kind=real_kind), dimension(: ,:), intent(out) :: edout, esout
                logical, intent(out) :: error
            end subroutine monrad
        end interface

        public :: real_kind, oasim_lib, calc_unit, pi, pi2, rad, rad_1
contains
        type(oasim_lib) function init_lib(filename, lat, lon, error) result(lib)
            implicit none

            character(len=*), intent(in) :: filename
            real(kind=real_kind), dimension(:), target, intent(in) :: lat, lon
            logical, intent(out) :: error

            type(table) :: original

            lib%init_parameters = parameters(filename, error)
            if (error) return

            original = table(lib%init_parameters%atmo_file, error, 6)
            if (error) return
            lib%atmo_adapted = table(lib%init_parameters%bin_file, error)
            if (error) goto 999
            call lib%atmo_adapted%from_table(original)
            call original%finalize()
            
            lib%rows = lib%atmo_adapted%get_rows()

            original = table(lib%init_parameters%abso_file, error, 2)
            if (error) return
            lib%abso_adapted = table(lib%init_parameters%bin_file, error)
            if (error) goto 999
            call lib%abso_adapted%from_table(original)
            call original%finalize()

            original = table(lib%init_parameters%slingo_file, error, 6)
            if (error) return
            lib%slingo_adapted = table(lib%init_parameters%bin_file, error)
            if (error) goto 999
            call lib%slingo_adapted%from_table(original)
            call original%finalize()

            error = .true.

            if (size(lat) /= size(lon)) return
            lib%m_size = size(lat)
            lib%lat => lat
            lib%lon => lon

            allocate(lib%up(lib%m_size, 3))
            allocate(lib%no(lib%m_size, 3))
            allocate(lib%ea(lib%m_size, 3))

            call lib%localvec() 

            error = .false.

            return

999         call original%finalize()
        end function init_lib

        subroutine finalize_lib(self)
            implicit none

            class(oasim_lib) :: self

            if (allocated(self%up)) deallocate(self%up)
            if (allocated(self%no)) deallocate(self%no)
            if (allocated(self%ea)) deallocate(self%ea)

            call self%atmo_adapted%finalize()
            call self%abso_adapted%finalize()
            call self%slingo_adapted%finalize()
        end subroutine finalize_lib

        subroutine localvec(self)
            implicit none
    
            class(oasim_lib) :: self

            real(kind=real_kind), dimension(:), allocatable :: cosy, upxy_1

            allocate(cosy(self%m_size))
            allocate(upxy_1(self%m_size))
    
            cosy = cos(self%lat * rad_1)
    
            if (self%init_parameters%local_time) then
                self%up(:,1) = cosy
                self%up(:,2) = 0.0
            else
                self%up(:,1) = cosy * cos(self%lon * rad_1)
                self%up(:,2) = cosy * sin(self%lon * rad_1)
            end if
            self%up(:,3) = sin(self%lat * rad_1)
    
            upxy_1 = 1.0d0 / sqrt(self%up(:,1) * self%up(:,1) &
                            + self%up(:,2) * self%up(:,2))
            self%ea(:,1) = - self%up(:,2) * upxy_1
            self%ea(:,2) = self%up(:,1) * upxy_1
            self%ea(:,3) = 0.0d0
            self%no(:,1) = - self%up(:,3) * self%ea(:,2)
            self%no(:,2) = self%up(:,3) * self%ea(:,1)
            self%no(:,3) = self%up(:,1) * self%ea(:,2) &
                            - self%up(:,2) * self%ea(:,1)

            if (allocated(cosy)) deallocate(cosy)
            if (allocated(upxy_1)) deallocate(upxy_1)
        end subroutine localvec

        type(calc_unit) function init_calc(p_size, lib) result(calc)
            implicit none

            integer, intent(in) :: p_size
            type(oasim_lib), target, intent(in) :: lib

            integer :: rows, i
            real(kind=real_kind) :: rlam, t, tlog, fac
            real(kind=real_kind), dimension(:), pointer :: aw, bw

            real(kind=real_kind), parameter :: a0 = 0.9976d0
            real(kind=real_kind), parameter :: a1 = 0.2194d0
            real(kind=real_kind), parameter :: a2 = 5.554d-2
            real(kind=real_kind), parameter :: a3 = 6.7d-3
            real(kind=real_kind), parameter :: b0 = 5.026d0
            real(kind=real_kind), parameter :: b1 = -1.138d-2
            real(kind=real_kind), parameter :: b2 = 9.552d-6
            real(kind=real_kind), parameter :: b3 = -2.698d-9            

            calc%lib => lib            
            rows = lib%rows
            calc%p_size = p_size

            aw => lib%abso_adapted%tab(:,1)
            bw => lib%abso_adapted%tab(:,2)

            allocate(calc%ed(rows))
            allocate(calc%es(rows))
            allocate(calc%td(rows))
            allocate(calc%ts(rows))
            allocate(calc%tcd(rows))
            allocate(calc%tcs(rows))
            allocate(calc%edclr(rows))
            allocate(calc%esclr(rows))
            allocate(calc%edcld(rows))
            allocate(calc%escld(rows))
            allocate(calc%tgas(rows))
            allocate(calc%ta(rows))
            allocate(calc%wa(rows))
            allocate(calc%asym(rows))
            allocate(calc%rlamu(rows))
            allocate(calc%wfac(rows))
            allocate(calc%ed_b(rows))
            allocate(calc%et_b(rows))
            allocate(calc%wv_b(rows))
            allocate(calc%rh_b(rows))

            allocate(calc%solz(p_size))
            allocate(calc%integ(p_size))
            allocate(calc%sunv(p_size))
            allocate(calc%sunn(p_size))
            allocate(calc%sune(p_size))
            allocate(calc%sunz(p_size))

            allocate(calc%rod(p_size, rows))
            allocate(calc%ros(p_size, rows))
            allocate(calc%eda(p_size, rows))
            allocate(calc%esa(p_size, rows))
            
            calc%tcd = 0.0
            calc%tcs = 0.0
            calc%td = 0.0
            calc%ts = 0.0
            do i = 1, rows
                calc%rlamu(i) = (lib%atmo_adapted%get_low(i) + lib%atmo_adapted%get_high(i)) * 0.5d-3
            end do

            do i = 1, rows
                rlam = (lib%atmo_adapted%get_low(i) + lib%atmo_adapted%get_high(i)) * 0.5
                if (rlam < 900.0) then
                    t = exp(-(aw(i) + 0.5 * bw(i)))
                    tlog = log(1.0d-36 + t)
                    fac = a0 + a1 * tlog + a2 * tlog * tlog + a3 * tlog * tlog * tlog
                    calc%wfac(i) = min(fac, 1.0d0)
                    calc%wfac(i) = max(fac, 0.0d0)
                else
                    fac = b0 + b1 * rlam + b2 * rlam * rlam + b3 * rlam * rlam * rlam
                    calc%wfac(i) = max(fac, 0.0d0)
                end if
            end do
        end function init_calc

        subroutine finalize_calc(self)
            implicit none

            class(calc_unit) :: self

            if (allocated(self%ed)) deallocate(self%ed)
            if (allocated(self%es)) deallocate(self%es)
            if (allocated(self%td)) deallocate(self%td)
            if (allocated(self%ts)) deallocate(self%ts)
            if (allocated(self%tcd)) deallocate(self%tcd)
            if (allocated(self%tcs)) deallocate(self%tcs)
            if (allocated(self%edclr)) deallocate(self%edclr)
            if (allocated(self%esclr)) deallocate(self%esclr)
            if (allocated(self%edcld)) deallocate(self%edcld)
            if (allocated(self%escld)) deallocate(self%escld)
            if (allocated(self%tgas)) deallocate(self%tgas)
            if (allocated(self%ta)) deallocate(self%ta)
            if (allocated(self%wa)) deallocate(self%wa)
            if (allocated(self%asym)) deallocate(self%asym)
            if (allocated(self%rlamu)) deallocate(self%rlamu)
            if (allocated(self%wfac)) deallocate(self%wfac)
            if (allocated(self%ed_b)) deallocate(self%ed_b)
            if (allocated(self%et_b)) deallocate(self%et_b)
            if (allocated(self%wv_b)) deallocate(self%wv_b)
            if (allocated(self%rh_b)) deallocate(self%rh_b)

            if (allocated(self%solz)) deallocate(self%solz)
            if (allocated(self%integ)) deallocate(self%integ)
            if (allocated(self%sunv)) deallocate(self%sunv)
            if (allocated(self%sunn)) deallocate(self%sunn)
            if (allocated(self%sune)) deallocate(self%sune)
            if (allocated(self%sunz)) deallocate(self%sunz)

            if (allocated(self%rod)) deallocate(self%rod)
            if (allocated(self%ros)) deallocate(self%ros)
            if (allocated(self%eda)) deallocate(self%eda)
            if (allocated(self%esa)) deallocate(self%esa)    
        end subroutine finalize_calc
end module oasim