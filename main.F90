program main

    use quad

    interface
    ! Prototype of the function implementation
    function fnmall(x,p)
        implicit none
        double precision     , dimension(:), intent(in) :: x,p
        double precision     , dimension(:) , allocatable :: fnmall
    end function fnmall

    end interface

    integer, dimension(:), allocatable :: grid
    integer  ,  dimension(:), allocatable :: n     
    double precision  ,  dimension(:), allocatable ::   r1, r2, p
    double precision :: tol, alpha, A, sigw,s
    double precision, dimension(3) :: beta,G,Na,sigma,Sm
    integer :: i, nm, np, d
    logical :: adap
    external f1, f2
    integer, parameter             :: nquad=16
    double precision, dimension(nquad) :: b,xi,wi
    double precision, dimension(2) :: endpts
    real :: pi
    pi = 4.*atan(1.)

    ! Input data as obtained from 3mode.txt
    ! alpha, A, G, beta and Sm are computed separately from the integrals
    alpha=5.537867855919976d-6
    A = 1.177280550774527d-7
    beta=(/2.982401212999331d0, 2.982401212999331d0, 2.982401212999331d0/)
    G=(/0.840528751639683,  0.840528751639683,  0.840528751639683/)
    G=G*1.0d-6
    Na=(/300.0d0, 200.0d0, 2.0d0/)
    sigma=(/1.8d0,   1.8d0, 1.8d0/)
    Sm=(/0.007773884286541d0, 0.001187233125088d0, 0.000037543607889d0/)
    sigw=20.0d0

    nm=3                     ! 3 modes
    allocate(n(1))
    allocate(grid(1))
    allocate(r1(1))
    allocate(r2(1))
    tol=1e-6                 ! tolerance (only needed in adaptive integration)
    grid=1                   ! Legendre quadrature
    n=9                      ! 9 quadrature points
    r1(1)=0.0d0              ! integration lower bound
    r2(1)=5.0d0*sigw         ! integration upper bound
    d = 1                    ! one dimensional integration
    np=5*nm+4                ! total number of parameters
    adap = .false.           ! no adaptivity

    ! Building input parameter vector p
    allocate(p(np))
    p(1)=alpha
    p(2)=A
    p(3:2+nm)=beta(1:nm)
    p(3+nm:2+2*nm)=G(1:nm)
    p(3+2*nm:2+3*nm)=Na(1:nm)
    p(3+3*nm:2+4*nm)=sigma(1:nm)
    p(3+4*nm:2+5*nm)=Sm(1:nm)
    p(np)=sigw

    write(*,*) '       Mode', '             fn   ', '                  fm   ',&
              '                 fluxfn   ', '                fluxfm   '

    ! quad_integral_PDF is called from the quadrature library
    do i=1,nm
      p(np-1)=i
       write(*,*) i, quad_integral_PDFv(fnmall,d,grid,p,r1,r2,adap,n,tol)
    enddo

    ! Test Hermite quadrature
    call gaussq(4,nquad,0.0d0,0.0d0,0,endpts,b,xi,wi)
    ! sum quadrature weights
    s = 0.0
    do i=1,nquad
        s = s + wi(i)
    enddo
    print *, "sum of Hermite weights = ", s/sqrt(pi)

    ! Test Legendre quadrature
    call gaussq(1,nquad,0.0d0,0.0d0,0,endpts,b,xi,wi)
    ! sum quadrature weights
    s = 0.0
    do i=1,nquad
        s = s + wi(i)
    enddo
    print *, "sum of Legendre weights = ", s/2.0

end program main

! Include the function implementation
#include "fnfm.F90"

!   Sample functions
    double precision function f1(x)
       double precision, intent(in) :: x
       f1 = cos(x)
    end function

    double precision function f2(x)
       double precision, intent(in) :: x
       f2 = (1 - tanh(x)**2)*dexp(x**2)
    end function
