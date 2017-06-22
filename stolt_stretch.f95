program stolt_stretch
! evaluates Stolt stretch parameter W for a set of depths for the given velocity array
! which can then be used as a migration parameter

! INPUT FILE: data.txt
!   first row: array length
!   then follows array of depths and array of velocities [m and m/s]
! 
!   by S. Vakeva, developed for the OpenFIRE group at Institute of Seismology,
!   University of Helsinki
!
!   Reference: Fomel, Vallant (2001): "Evaluating the Stolt-stretch parameter"
!              Journal of Seismic Exploration, v. 9, 319-335

    implicit none

    real, dimension(:), allocatable :: z,vz,vrms,tv,s,ss, &
                                        W,vint,temp,eta
    real :: v0
    integer :: n,i

    open(unit=15, file='data.txt', status='old', &
        action='read')
    read (15,*) n
    allocate(z(n),vz(n),vrms(n),vint(n-1),tv(n),s(n),ss(n), &
        W(n),temp(n),eta(n))
    read(15,*) z
    read(15,*) vz
    close(15)

    ! z to meters
    z=z*1000.0

    ! calculate twt
    tv(1)=0.0
    do i=2,n
        tv(i)=tv(i-1)+(z(i)-z(i-1))/(vz(i)-vz(i-1))*log(vz(i)/vz(i-1))
    end do
    ! calculate vint
    do  i=1,n-1
        vint(i)=(z(i+1)-z(i))/(tv(i+1)-tv(i))
    end do
    ! calculate vrms
    vrms(1)=vz(1)**2.0
    do i=2,n
        vrms(i)=(tv(i)-tv(i-1))*vint(i-1)**2.0+tv(i-1)*vrms(i-1)
        vrms(i)=vrms(i)/tv(i)
    end do

    ! calculate parameter of heterogeneity S
    ss(1)=1.0
    temp(1)=0.0
    do i=2,n
        temp(i)=temp(i-1)+(vz(i)**5.0-vz(i-1)**5.0)*(tv(i)-tv(i-1)) &
            /5.0/(vz(i)-vz(i-1))
        ss(i)=temp(i)/vrms(i)**2.0/tv(i)
    end do

    ! calculate eta
    eta(1)=0.0
    do i=2,n
        eta(i)=eta(i-1)+0.5*(vz(i)+vz(i-1))*(z(i)-z(i-1))
    end do

    ! calculate stolt stretch
    s(1)=0.0
    temp(1)=0.0
    do i=2,n
        temp(i)=temp(i-1)+0.5*(eta(i)+eta(i-1))*(tv(i)-tv(i-1))
        s(i)=sqrt(2.0*temp(i))
    end do

    ! calculate W
    W(1)=1.0
    do i=2,n
        W(i)=1-s(i)**2.0/vrms(i)/tv(i)**2.0*(vz(i)**2.0/vrms(i)-ss(i))
    end do


    ! print velocities
    write(*,*) 'Velocity array:'
    write(*,100) ' z/km  ',' twt/s ','  v_z  ',' vint  ',' vrms  ','   W   '
    100 format(5a7,3x,a6)
    do i=1,n
        write(*,101) z(i)/1000.0,2.0*tv(i),vz(i),vint(i),sqrt(vrms(i)),W(i)
    end do
    101 format(5f7.1,3x,f6.4)

end program stolt_stretch
