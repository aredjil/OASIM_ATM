module oasim_common
    use :: iso_c_binding

    implicit none
    
    integer, parameter :: string_length = 1024
    integer, parameter :: real_kind = c_double
    type:: light_return 
        real(kind=real_kind), allocatable:: tcd(:)
        real(kind=real_kind), allocatable:: tcs(:)
    end type light_return
end module oasim_common
