program benchmark
    use oasim
    ! Calling the netcdf library to use it for opening and reading the data from the file 
    use netcdf 
    use openacc
    
    implicit none
    !----------------------------
    logical :: error
    ! Vriables for timing the main loop 
    !----------------------------
     real :: start, finish, elapsed
    ! netCDF variables section
    !------------------------------------------------------------------------------------ 
    ! Netcdf variables 
    !-----------------
    integer:: ncid, varid, retval
    ! Dimensions variables  412300
    !---------------------
    integer, parameter :: n_iter = 300, n_points = 1000, n_waves = 33

    
    ! Variables to read from NetCDF
    !-------------------
    integer, allocatable :: iyr(:), iday(:)
    real(kind=real_kind), allocatable ::  sec_b(:), sec_e(:)
    real(kind=real_kind), allocatable ::  sp(:,:), msl(:,:), ws10(:,:), tco3(:,:), t2m(:,:)
    real(kind=real_kind), allocatable ::  d2m(:,:), tcc(:,:), tclw(:,:), cdrem(:,:)
    real(kind=real_kind), allocatable ::  taua(:,:,:), asymp(:,:,:), ssalb(:,:,:)
    real(kind=real_kind), allocatable ::  lat(:), lon(:)

    ! Variables for OASIM computation
    !--------------------------------
    integer, dimension(:), allocatable :: points
    real(kind=real_kind), dimension(:), allocatable :: slp, wsm, oz, wv, rh, ccov, rlwp, cdre
    real(kind=real_kind), dimension(:,:), allocatable :: edout, esout
    real(kind=real_kind), dimension(:,:), allocatable :: taua_slice, asymp_slice, ssalb_slice
    
    ! Physical constants
    !-------------------
    real(kind=real_kind), parameter :: a1 = 611.21d0 ! Pascal
    real(kind=real_kind), parameter :: a3 = 17.502d0 ! dimensionless
    real(kind=real_kind), parameter :: a4 = 32.19d0  ! Kelvin
    real(kind=real_kind), parameter :: To = 273.16d0 ! Kelvin
    real(kind=real_kind), parameter :: b1 = 0.14d0 * 0.01d0  ! cm/Pascal
    real(kind=real_kind), parameter :: b2 = 0.21d0   ! cm

    ! Working variables
    !------------------
    real(kind=real_kind) :: T, Td, es_Td, es_T
    integer :: i, j, iter
    ! File path 
    !----------
    character(len=256) :: home, FILE_NAME
    ! OASIM objects
    !--------------
    type(oasim_lib) :: lib
    type(calc_unit) :: calc

    ! Allocating memory for the variables 
    !------------------------------------
    allocate(iyr(n_iter), iday(n_iter))
    allocate(sec_b(n_iter), sec_e(n_iter))
    allocate(sp(n_points, n_iter), msl(n_points, n_iter), ws10(n_points, n_iter))
    allocate(tco3(n_points, n_iter), t2m(n_points, n_iter), d2m(n_points, n_iter))
    allocate(tcc(n_points, n_iter), tclw(n_points, n_iter), cdrem(n_points, n_iter))
    allocate(taua(n_points, n_waves, n_iter))
    allocate(asymp(n_points, n_waves, n_iter))
    allocate(ssalb(n_points, n_waves, n_iter))
    allocate(lat(n_points), lon(n_points))


    ! Opening the netcdf file and reading data
    !-----------------------------------------
    call get_environment_variable("HOME", home)
    FILE_NAME = trim(home)//"/projects/OASIM_ATM/data000.nc"
    call check(nf90_open(FILE_NAME, NF90_NOWRITE, ncid))
    
    ! Reading time variables
    call check(nf90_inq_varid(ncid, "iyr", varid))
    call check(nf90_get_var(ncid, varid, iyr))
    call check(nf90_inq_varid(ncid, "iday", varid))
    call check(nf90_get_var(ncid, varid, iday))
    ! Reading the start and end seconds of the day     
    ! sec_b = 28800.0d0  ! 8 AM in seconds
    ! sec_e = 36000.0d0  ! 10 AM in seconds
    call check(nf90_inq_varid(ncid, "sec_b", varid))
    call check(nf90_get_var(ncid, varid, sec_b))
    call check(nf90_inq_varid(ncid, "sec_e", varid))
    call check(nf90_get_var(ncid, varid, sec_e))


    ! call check(nf90_inq_varid(ncid, "lat", varid))
    ! call check(nf90_get_var(ncid, varid, lat))
    ! call check(nf90_inq_varid(ncid, "lon", varid))
    ! call check(nf90_get_var(ncid, varid, lon))

    ! Reading atmospheric variables
    call check(nf90_inq_varid(ncid, "sp", varid))
    call check(nf90_get_var(ncid, varid, sp))
    call check(nf90_inq_varid(ncid, "msl", varid))
    call check(nf90_get_var(ncid, varid, msl))
    call check(nf90_inq_varid(ncid, "ws10", varid))
    call check(nf90_get_var(ncid, varid, ws10))
    call check(nf90_inq_varid(ncid, "tco3", varid))
    call check(nf90_get_var(ncid, varid, tco3))
    call check(nf90_inq_varid(ncid, "t2m", varid))
    call check(nf90_get_var(ncid, varid, t2m))
    call check(nf90_inq_varid(ncid, "d2m", varid))
    call check(nf90_get_var(ncid, varid, d2m))
    call check(nf90_inq_varid(ncid, "tcc", varid))
    call check(nf90_get_var(ncid, varid, tcc))
    call check(nf90_inq_varid(ncid, "tclw", varid))
    call check(nf90_get_var(ncid, varid, tclw))
    call check(nf90_inq_varid(ncid, "cdrem", varid))
    call check(nf90_get_var(ncid, varid, cdrem))
    call check(nf90_inq_varid(ncid, "taua", varid))
    call check(nf90_get_var(ncid, varid, taua))
    call check(nf90_inq_varid(ncid, "asymp", varid))
    call check(nf90_get_var(ncid, varid, asymp))
    call check(nf90_inq_varid(ncid, "ssalb", varid))
    call check(nf90_get_var(ncid, varid, ssalb))

    ! Close the NetCDF file
    call check(nf90_close(ncid))

    ! write(*, *) "NetCDF data read successfully"
    ! write(*, *) "Number of points: ", n_points
    ! write(*, *) "Number of iterations: ", n_iter
    ! write(*, *) "Number of wavelengths: ", n_waves

    ! Initialize OASIM library
    !-------------------------
    lib = oasim_lib("config.yaml", lat, lon, error)
    
    if (error) then
        write(*, *) "Error initializing OASIM library"
        stop 1
    end if

    ! write(*, *) "OASIM library initialized successfully"

    ! Allocate working arrays
    !------------------------
    allocate(points(n_points))
    allocate(slp(n_points), wsm(n_points), oz(n_points))
    allocate(wv(n_points), rh(n_points), ccov(n_points))
    allocate(rlwp(n_points), cdre(n_points))
    allocate(edout(n_points, lib%rows), esout(n_points, lib%rows))
    allocate(taua_slice(n_points, lib%rows))
    allocate(asymp_slice(n_points, lib%rows))
    allocate(ssalb_slice(n_points, lib%rows))

    ! Initialize points array
    do i = 1, n_points
        points(i) = i
    end do

    ! Initialize calculation unit
    calc = calc_unit(n_points, lib)

    ! Main computation loop
    !---------------------
    
     call cpu_time(start) ! Starting the timing 
    ! write(*, *) "Starting benchmark computations..."
    
    main_loop:do iter = 1, n_iter
        ! write(*, '(A,I0,A,I0)') "Processing iteration ", iter, " of ", n_iter
        
        ! Convert units and calculate derived variables for all points
        do i = 1, n_points
            ! Unit conversions
            slp(i) = sp(i, iter) / 100.0d0    ! Pa to hPa
            wsm(i) = ws10(i, iter)            ! m/s
            oz(i) = tco3(i, iter)             ! Total column ozone
            
            ! Calculate relative humidity and precipitable water
            T = t2m(i, iter)
            Td = d2m(i, iter)
            es_Td = a1 * exp(a3 * (Td - To) / (Td - a4))
            es_T = a1 * exp(a3 * (T - To) / (T - a4))
            rh(i) = 100.0d0 * es_Td / es_T
            wv(i) = b1 * es_Td * (sp(i, iter)/100.0d0) / (msl(i, iter)/100.0d0) + b2
            
            ! Cloud properties
            ccov(i) = tcc(i, iter)
            rlwp(i) = tclw(i, iter) * 1000.0d0  ! Convert to g/m²
            cdre(i) = cdrem(i, iter)
            
            ! Extract optical properties for this iteration
            taua_slice(i, :) = taua(i, :, iter)
            asymp_slice(i, :) = asymp(i, :, iter)
            ssalb_slice(i, :) = ssalb(i, :, iter)
        end do
        
        ! Call OASIM monrad function
        call calc%monrad(points, iyr(iter), iday(iter), sec_b(iter), sec_e(iter), &
                        sp(:, iter), msl(:, iter), ws10(:, iter), tco3(:, iter), &
                        t2m(:, iter), d2m(:, iter), tcc(:, iter), tclw(:, iter), &
                        cdrem(:, iter), taua_slice, asymp_slice, ssalb_slice, &
                        edout, esout, error)
        if (error) then
            write(*, *) "Error in OASIM calculation at iteration ", iter
            stop 1
        end if
        
    end do main_loop
    ! End the timing 
     call cpu_time(finish)

     elapsed = finish - start
     print *, 'Elapsed CPU time (seconds):', elapsed
    ! Clean up
    !---------
    call calc%finalize()
    call lib%finalize()
    
    deallocate(iyr, iday)
    deallocate(sec_b, sec_e)
    deallocate(sp, msl, ws10, tco3, t2m)
    deallocate(d2m, tcc, tclw, cdrem)
    deallocate(taua, asymp, ssalb)
    deallocate(lat, lon)
    deallocate(points)
    deallocate(slp, wsm, oz, wv, rh, ccov, rlwp, cdre)
    deallocate(edout, esout)
    deallocate(taua_slice, asymp_slice, ssalb_slice)

    contains
        ! Simple subroutine to check the status of reading the netcdf file 
        subroutine check(status)
            integer, intent(in) :: status
            
            if(status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
        end subroutine check  
end program benchmark
