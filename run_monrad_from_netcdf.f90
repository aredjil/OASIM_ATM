program run_monrad_from_netcdf
    use netcdf
    use oasim, only: calc_unit, real_kind
    use oasim_monrad
    implicit none

    type(calc_unit) :: cu
    integer :: ncid, varid, retval
    integer :: npoints, nrows
    integer, allocatable :: points(:)
    integer :: iyr, iday
    real(kind=real_kind) :: sec_b, sec_e
    real(kind=real_kind), allocatable :: sp(:), msl(:), ws10(:), tco3(:), t2m(:), d2m(:), tcc(:), tclw(:), cdrem(:)
    real(kind=real_kind), allocatable :: taua(:,:), asymp(:,:), ssalb(:,:)
    real(kind=real_kind), allocatable :: edout(:,:), esout(:,:)
    logical :: error

    ! Open the NetCDF file
    retval = nf90_open("data000.nc", NF90_NOWRITE, ncid)
    if (retval /= nf90_noerr) stop "Could not open data000.nc"

    ! Read dimensions
    call get_dim(ncid, "points", npoints)
    call get_dim(ncid, "rows", nrows)

    ! Allocate arrays
    allocate(points(npoints))
    allocate(sp(npoints), msl(npoints), ws10(npoints), tco3(npoints), t2m(npoints), d2m(npoints), tcc(npoints), tclw(npoints), cdrem(npoints))
    allocate(taua(npoints, nrows), asymp(npoints, nrows), ssalb(npoints, nrows))
    allocate(edout(npoints, nrows), esout(npoints, nrows))

    ! Read variables
    call get_var(ncid, "points", points)
    call get_var(ncid, "iyr", iyr)
    call get_var(ncid, "iday", iday)
    call get_var(ncid, "sec_b", sec_b)
    call get_var(ncid, "sec_e", sec_e)
    call get_var(ncid, "sp", sp)
    call get_var(ncid, "msl", msl)
    call get_var(ncid, "ws10", ws10)
    call get_var(ncid, "tco3", tco3)
    call get_var(ncid, "t2m", t2m)
    call get_var(ncid, "d2m", d2m)
    call get_var(ncid, "tcc", tcc)
    call get_var(ncid, "tclw", tclw)
    call get_var(ncid, "cdrem", cdrem)
    call get_var(ncid, "taua", taua)
    call get_var(ncid, "asymp", asymp)
    call get_var(ncid, "ssalb", ssalb)

    ! Set up calc_unit object (you may need to initialize more fields)
    cu%p_size = npoints
    cu%lib%rows = nrows

    ! Call monrad
    call monrad(cu, points, iyr, iday, sec_b, sec_e, sp, msl, ws10, tco3, t2m, d2m, &
                tcc, tclw, cdrem, taua, asymp, ssalb, edout, esout, error)

    if (error) then
        print *, "monrad returned error"
    else
        print *, "monrad executed successfully"
        ! Optionally, print or process edout, esout
    end if

    call nf90_close(ncid)

contains

    subroutine get_dim(ncid, name, dimlen)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: dimlen
        integer :: dimid, retval
        retval = nf90_inq_dimid(ncid, name, dimid)
        if (retval /= nf90_noerr) stop "Dimension "//trim(name)//" not found"
        retval = nf90_inq_dimlen(ncid, dimid, dimlen)
        if (retval /= nf90_noerr) stop "Could not get length of dimension "//trim(name)
    end subroutine get_dim

    subroutine get_var(ncid, name, var)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        ! Generic interface for different types
        ! Integer scalar
        integer, optional, intent(out) :: var
        ! Real(kind=real_kind) scalar
        real(kind=real_kind), optional, intent(out) :: var
        ! Integer array
        integer, optional, intent(out), dimension(:) :: var
        ! Real(kind=real_kind) array
        real(kind=real_kind), optional, intent(out), dimension(:) :: var
        ! Real(kind=real_kind) 2D array
        real(kind=real_kind), optional, intent(out), dimension(:,:) :: var
        integer :: varid, retval
        if (present(var) .and. size(shape(var)) == 0) then
            retval = nf90_inq_varid(ncid, name, varid)
            if (retval /= nf90_noerr) stop "Variable "//trim(name)//" not found"
            retval = nf90_get_var(ncid, varid, var)
            if (retval /= nf90_noerr) stop "Could not read variable "//trim(name)
        else if (present(var) .and. size(shape(var)) == 1) then
            retval = nf90_inq_varid(ncid, name, varid)
            if (retval /= nf90_noerr) stop "Variable "//trim(name)//" not found"
            retval = nf90_get_var(ncid, varid, var)
            if (retval /= nf90_noerr) stop "Could not read variable "//trim(name)
        else if (present(var) .and. size(shape(var)) == 2) then
            retval = nf90_inq_varid(ncid, name, varid)
            if (retval /= nf90_noerr) stop "Variable "//trim(name)//" not found"
            retval = nf90_get_var(ncid, varid, var)
            if (retval /= nf90_noerr) stop "Could not read variable "//trim(name)
        end if
    end subroutine get_var

end program run_monrad_from_netcdf