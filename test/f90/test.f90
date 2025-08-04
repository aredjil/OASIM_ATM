program test
    use oasim
    implicit none

    logical :: error
    integer :: year, day !,i
    real(kind=real_kind) :: sec_b, sec_e
    integer, parameter :: m_size = 1
    integer, parameter :: p_size = 1
    real(kind=real_kind), parameter :: a1 = 611.21d0 ! Pascal
    real(kind=real_kind), parameter :: a3 = 17.502d0 ! dimensionless
    real(kind=real_kind), parameter :: a4 = 32.19d0  ! Kelvin
    real(kind=real_kind), parameter :: To = 273.16d0 ! Kelvin
    real(kind=real_kind), parameter :: b1 = 0.14d0 * 0.01d0  ! cm/Pascal
    real(kind=real_kind), parameter :: b2 = 0.21d0   ! cm
    real(kind=real_kind)                            :: sp, msl,d2m,t2m,tcc,ws10,tclw,tco3,cdrem,cldtcm,T,Td,es_Td,es_T
    real(kind=real_kind), dimension(:), allocatable :: lat, lon
    integer, dimension(:), allocatable :: points
    real(kind=real_kind), dimension(:), allocatable :: slp, wsm, oz, wv, rh, ccov, rlwp, cdre
    real(kind=real_kind), dimension(:,:), allocatable :: taua, asymp, ssalb, edout, esout
    type(oasim_lib) :: lib
    type(calc_unit) :: calc

    allocate(lat(m_size))
    allocate(lon(m_size))

    ! do i = 1, m_size
    !     lat(i) = float(i) * 0.25d0 - 3.0d0
    !     lon(i) = 40.0d0 + float(i) * 0.10d0
    ! end do

    lon=7.2916665
    lat=42.9375

    lib = oasim_lib("config.yaml", lat, lon, error)

    if (error) then
        write(*, *) "Something wrong happened."
        stop 1
    end if

    write(*, *) "Configuration file correctly read:"
    write(*, *) " "
    write(*, *) trim(lib%init_parameters%atmo_file)
    write(*, *) trim(lib%init_parameters%abso_file)
    write(*, *) trim(lib%init_parameters%slingo_file)
    write(*, *) lib%init_parameters%integration_step_secs
    write(*, *) lib%init_parameters%max_integration_steps
    write(*, *) lib%init_parameters%zenith_avg
    write(*, *) lib%init_parameters%local_time
    write(*, *) trim(lib%init_parameters%bin_file)
    write(*, *) " "
    write(*, *) " "
    write(*, *) " "
    write(*, *) "Tables read correctly."
    write(*, *) " "
    write(*, *) " "
    write(*, *) "ATMO"
    write(*, *) " "
    call lib%atmo_adapted%print()
    write(*, *) " "
    write(*, *) " "
    write(*, *) "ABSO"
    write(*, *) " "
    call lib%abso_adapted%print()
    write(*, *) " "
    write(*, *) " "
    write(*, *) "SLINGO"
    write(*, *) " "
    call lib%slingo_adapted%print()
    write(*, *) " "
    write(*, *) " "
    ! write(*, *) lib%m_size
    ! write(*, *) " "
    ! do i = 1, lib%m_size
    !     write(*, *) lib%lat(i), lib%lon(i), lib%up(i,:), lib%no(i,:), lib%ea(i,:)
    ! end do

    allocate(slp(p_size))
    allocate(wsm(p_size))
    allocate(oz(p_size))
    allocate(wv(p_size))
    allocate(rh(p_size))
    allocate(ccov(p_size))
    allocate(rlwp(p_size))
    allocate(cdre(p_size))
    allocate(points(p_size))

    ! write(*, *) lib%rows, p_size

    allocate(taua(p_size, lib%rows))
    allocate(asymp(p_size, lib%rows))
    allocate(ssalb(p_size, lib%rows))
    allocate(edout(p_size, lib%rows))
    allocate(esout(p_size, lib%rows))

    points = 1

    sp=102377.13
    msl=102410.06
    d2m=284.10672
    t2m=286.96487
    tcc=12.03943
    ws10=10.0
    tclw=0.05
    tco3=0.0065
    cdrem=10.042508
    cldtcm=10.042508
    
    taua(1,:)=[0.12861817,0.11693237,0.11303711,0.10914185,0.10524658,0.10135132 &
    ,0.09745605,0.0934308,0.08906767,0.08519087,0.08188277,0.079188 &
    ,0.07697173,0.07506739,0.07332103,0.07159076,0.06982587,0.06804233 &
    ,0.0645105,0.05960191,0.05457829,0.05089346,0.0479044,0.04534442 &
    ,0.04311129,0.04106918,0.03914813,0.0373441,0.03567282,0.03332167 &
    ,0.02901247,0.01838189,0.00623265]
    asymp(1,:)=[0.78112185,0.76577055,0.76065344,0.7555364,0.75041926,0.74530214 &
    ,0.7401851,0.7352586,0.7308154,0.72576404,0.71973574,0.7124992 &
    ,0.7045832,0.69671303,0.68963724,0.68404067,0.6799671,0.6770385 &
    ,0.6733678,0.6695107,0.65929943,0.64360243,0.62599105,0.6091321 &
    ,0.59477234,0.58193725,0.5697906,0.5578721,0.54593724,0.5280586 &
    ,0.49310532,0.41066745,0.3164527]
    ssalb(1,:)=[0.9653854,0.9685926,0.96966165,0.9707307,0.9717998,0.97286886 &
    ,0.9739379,0.9750599,0.97628933,0.97734255,0.97821194,0.9788913 &
    ,0.9793933,0.97973263,0.9799186,0.9799566,0.9798687,0.9796886 &
    ,0.9791781,0.9784721,0.9785546,0.9787444,0.97678363,0.96918947 &
    ,0.95303714,0.9326887,0.91252154,0.897008,0.88944095,0.88677096 &
    ,0.89569014,0.89257544,0.88901585]

! Conversions
    slp = sp /100.0D0   
    wsm = ws10        ! m/s
    oz  = tco3
    
 ! relative humidity and
 ! precipitable water (absorption by water vapour)  
 
    
    T            = t2m
    Td           = d2m 
    es_Td        = a1 * exp( a3 * (Td-To)/(Td-a4) )
    es_T         = a1 * exp( a3 *  (T-To)/(T-a4) )
    rh           = 100. * es_Td /es_T
    wv           = b1 * es_Td * (sp/100.)/(msl/100.) + b2
    
    ccov = tcc
    
    rlwp = tclw*1000.
    
    cdre = cdrem
    
   ! cldtc  = cldtcm
    
    year = 2019
    day = 1
    sec_b = 28800.0d0
    sec_e = 36000.0d0
    calc = calc_unit(p_size, lib)
    call calc%monrad_debug(points, year, day, sec_b, sec_e, slp, wsm, oz, wv, & 
                          rh, ccov, rlwp, cdre, taua, asymp, ssalb, edout, esout, error)
    ! write(*, *) error
    write(*, *) "FINAL OUTPUT - direct"
    write(*, *) edout
    write(*, *) "FINAL OUTPUT - diffuse"
    write(*, *) esout      
                     
    call calc%monrad(points, year, day, sec_b, sec_e, [sp], [msl], [ws10], [tco3], [t2m], [d2m], &
                     [tcc], [tclw], [cdrem], taua, asymp, ssalb, edout, esout, error)

    ! write(*, *) error
    write(*, *) "FINAL OUTPUT - direct"
    write(*, *) edout
    write(*, *) "FINAL OUTPUT - diffuse"
    write(*, *) esout  


    call calc%finalize()
    deallocate(slp)
    deallocate(wsm)
    deallocate(oz)
    deallocate(wv)
    deallocate(rh)
    deallocate(ccov)
    deallocate(rlwp)
    deallocate(cdre)
    deallocate(points)
    deallocate(taua)
    deallocate(asymp)
    deallocate(ssalb)
    deallocate(edout)
    deallocate(esout)
    call lib%finalize()
    if (allocated(lat)) deallocate(lat)
    if (allocated(lon)) deallocate(lon)
end program test