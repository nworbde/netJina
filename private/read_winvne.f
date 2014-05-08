program read_winvne
    use iso_fortran_env, only: iostat_end
    implicit none

    character(len=*),parameter :: filename = 'winvne_v2.0.dat'
    character(len=*), parameter :: pfcn_fmt = '(23f3.2,f3.1)'
    character(len=*),parameter :: iso_fmt = '(a5)'
    character(len=*),parameter :: line0_fmt = '(a5,f12.3,2i4,f6.1,f10.3,a6)'
    character(len=5) :: iso_name
    character(len=6) :: provenance
    real :: A, S, Q
    integer :: Z, N
    integer :: i, iso_count,ios
    real, dimension(24) :: pfcn_temps, pfcn
    
    open(unit=10,file=filename,status='old',action='read')
    read(10,pfcn_fmt) pfcn_temps(:)
    write(*,'(f6.3)') pfcn_temps(:)

    A = 0.0
    S = 0.0
    Q = 0.0
    Z = 0
    N = 0
    iso_count = 0
    
    do
        ! try to read a full data line; if that fails, then read in the isotope 
        ! name
        read(10,fmt=line0_fmt,iostat=ios) iso_name,A,Z,N,S,Q,provenance
        if (ios /= 0 .or. A == 0.0) then
            backspace 10
            read(10,fmt=iso_fmt,iostat=ios) iso_name
            if (ios /= 0) exit
            iso_count = iso_count+1
            print *,iso_count,iso_name
        else
            read(10,*) pfcn(1:8)
            read(10,*) pfcn(9:16)
            read(10,*) pfcn(17:24)
        end if
    end do
    print *, 'obtained ',iso_count,' isotopes'
    close(10)

end program read_winvne
