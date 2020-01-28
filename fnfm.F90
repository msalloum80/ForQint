function fnmall(x,p)
    double precision, dimension(:), intent(in) :: x
    double precision, dimension(:), intent(in) :: p
    double precision     , dimension(:) ,allocatable :: fnmall
    integer :: np, nm, iout, i
    double precision :: alpha,A,xx,Smax,sigw,w,eta,zeta
    double precision  ,  dimension(:), allocatable :: beta,G,Na,sigma,Sm

    w=x(1)
    np=size(p)
    nm=(np-4)/5
    iout=p(np-1)
    if ((iout.gt.nm).or.(iout.lt.0)) then
     write(*,*) 'invalid paramters input or mode index'
     stop
    end if

    ! Inflating the paramters vector p
    alpha=p(1)
    A=p(2)
    allocate(beta(nm))
    allocate(G(nm))
    allocate(Na(nm))
    allocate(sigma(nm))
    allocate(Sm(nm))
    beta(1:nm)=p(3:2+nm)
    G(1:nm)=p(3+nm:2+2*nm)
    Na(1:nm)=p(3+2*nm:2+3*nm)
    sigma(1:nm)=p(3+3*nm:2+4*nm)
    Sm(1:nm)=p(3+4*nm:2+5*nm)
    sigw=p(np)

    ! Computing Smax
    Smax=0.0d0
    do i=1,nm    
     eta=2.0d0*(alpha*w)**1.5d0/Na(i)/dsqrt(G(i))/beta(i)
     zeta=2.0d0*A*dsqrt(alpha*w)/3.0d0/dsqrt(G(i))
     Smax=Smax+(0.5d0*dexp(2.5d0*dlog(sigma(i))**2)/Sm(i)**2*(zeta/eta)**1.5d0+ &
        (1.0d0+0.25d0*dlog(sigma(i)))*(Sm(i)**2/(eta+3*zeta))**0.75d0/Sm(i)**2)
    enddo
    Smax=1.0d0/dsqrt(Smax)

    xx=dlog((Sm(iout)/Smax)**2)/dlog(sigma(iout))/3.0d0/dsqrt(2.0d0)

    ! Building function output vector  
    allocate(fnmall(4))

    ! fn
    fnmall(1)=(1.0d0-derf(xx))*dexp(-(w/sigw)**2/2.0d0)/ &
      sqrt(2.0d0*3.1415926535897932d0)/sigw

    ! fm
    fnmall(2)=(1.0d0-derf(xx-1.4d0*dsqrt(2.0d0)*dlog(sigma(iout))))* &
       dexp(-(w/sigw)**2/2.0d0)/sqrt(2.0d0*3.1415926535897932d0)/sigw

    ! fluxfn
    fnmall(3)=fnmall(1)*0.5d0*Na(iout)*w

    ! fluxfm
    fnmall(4)=fnmall(2)*0.5d0*w
end function fnmall
