module quad
use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float,c_f_pointer,c_null_char, c_long_double, c_loc

implicit none

interface

! Function prototype with input parameters as a vector
double precision function myfunc(x,p)
    implicit none
    double precision     , dimension(:), intent(in) :: &
        x, & ! vector of input variables
        p    ! p = vector of parameters
end function

! Function prototype with input parameters as a vector and vector output
function myfuncd(x,p)
    implicit none
    double precision     , dimension(:), intent(in) :: &
        x, & ! vector of input variables
        p    ! p = vector of parameters
    double precision     , dimension(:) ,allocatable :: myfuncd
end function


end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! MAIN QUAD FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine quadcreate_aniso_full(grid,n,x,w,table,xb,wb)
    ! grid = The quadrature rule in each dimension (size = d)
    ! n = The number of quadrature points in each dimension (size = d)
    ! w = Vector of quadrature weigths of (size = Np)
    ! x = Matrix of quadrature coordinates of (size = d x Np)
    ! table = Matrix of quadrature table of (size = d x Np)
    ! wb = 3D Array of quadrature weigths of (size = Np x Np x 2)
    ! xb = 3D Array of quadrature points of (size = Np x Np x 2)
    integer, dimension(:), allocatable, intent(in) :: grid  
    integer  ,  dimension(:), allocatable, intent(in) :: n     
    double precision  ,  dimension(:), allocatable, intent(out) :: w    
    double precision  ,  dimension(:,:), allocatable, intent(out) :: x    
    integer  ,  dimension(:,:), allocatable, intent(out) :: table
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: wb    
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: xb
    integer :: nm,d,i,j,cnt,ntot
    double precision     ,  dimension(:,:), allocatable ::  ww,xx,xt
    integer  ,  dimension(:), allocatable ::  ind
    integer :: fs = 0
    integer(c_int)                                   :: Quadd
    double precision :: pi,cx,cw
    integer, dimension(:), allocatable :: n1
    integer, dimension(:), allocatable :: r1
    double precision, dimension(:), allocatable :: work,xi,wi 
    double precision, dimension(2) :: endpts ! use standard endpoints

    allocate(n1(1))
    allocate(r1(1))
    pi = 4.0d0*datan(1.0d0)
    cx = sqrt(2.0d0)
    cw = sqrt(pi)

    d=size(n,1)
    if (d.ne.size(grid,1)) then
        write(0,*) 'Error: the size of n and grid_type should be the same as ndim'
        stop
    end if
      
    nm = maxval(n,n.gt.1)
    allocate(xx(d,nm))
    allocate(ww(d,nm))

    if (present(xb).and.present(wb)) then
        do i=1,d
            allocate(xt(1,n(i)))
            ww(i,1:n(i)) = wb(n(i),1:n(i),grid(i))
            xt(1,:) = xb(n(i),1:n(i),grid(i))    
            xx(i,1:n(i)) = xt(1,:)
            deallocate(xt)
        enddo
    else
    ! Creates a 1D quadrature rule; extract points and weights from it,
    ! and use them to build a multi-D quadrature rule, one dimension at the time.
    ! (Could also do multi-D quadrature directly.)
    ! This would be a good place to call gaussq instead of Trilinos-based routines.
    r1(1)=grid(1)
    n1(1)=n(1)

    allocate(xi(n1(1)))
    allocate(wi(n1(1)))
    allocate(work(n1(1)))
    ! call gaussq(r1(1),n1(1),0.0d0,0.0d0,0,endpts,work,xi,wi)
    do i=1,d
        r1(1)=grid(i)
        n1(1)=n(i)
        if (r1(1) == 1) then
            call gaussq(1,n1(1),0.0d0,0.0d0,0,endpts,work,xi,wi)
            ww(i,1:n(i)) = .5*wi ! Normalized Gauss-Legendre
        end if
        if (r1(1) == 2) then
            call gaussq(4,n1(1),0.0d0,0.0d0,0,endpts,work,xi,wi)
            ww(i,1:n(i)) = wi/sqrt(pi) ! normalized Gauss-Hermite
        end if
        xx(i,1:n(i)) = xi
    enddo
        deallocate(xi)
        deallocate(wi)
        deallocate(work)
    end if

    ntot=1
    do i=1,d
        ntot=ntot*n(i)
    enddo
    allocate(ind(d))
    do i=1,d
        ind(i)=1
    enddo
    cnt=1
    allocate(x(d,ntot))
    allocate(table(d,ntot))
    allocate(w(ntot))
    w = 1.0d0
    do while (cnt.le.ntot)
        do j=1,d
            x(j,cnt) = xx(j,ind(j))
            w(cnt) = w(cnt)*ww(j,ind(j))
            table(j,cnt) = ind(j)
        enddo
        ind(1)=ind(1)+1
        i = 1
        do while ((i.le.(d-1)).and.((ind(i)-1).eq.n(i)))
            ind(i)=1
            i=i+1
            ind(i)=ind(i)+1
        enddo
        cnt=cnt+1
    end do
end subroutine quadcreate_aniso_full

! Function that integrates a function f according to a given distribution
! The function to be integrated has input parameters p
double precision function quad_integral(func,d,q,p,r1,r2,adap,Ni,tol,t,pmu,psigma,xb,wb)
    ! func = function to be integrated
    ! d = number of dimensions
    ! t = vector of strings denoting the type of probability distributions in each dimension (size = d)
    ! q = vector of strings denoting the type of quadrature in each dimension (size = d)
    ! p = vector of function parameters
    ! pmu = vector of distribution means (size = d)
    ! psigma = cholesky decomposition of the variance matrix (size = dxd)
    ! r1,r2 = vectors of the function variables ranges (size = d) [only used in Legendre quandratures]
    ! adap = integer denoting if adaptivity is enabled (adap=0 means no adaptivity)
    ! tol = relative tolerance [only used when adaptivity is enabled]
    ! wb = 3D Array of quadrature weigths of (size = n x n x 2) (optional, adding wb and xb accelerates the integration)
    ! xb = 3D Array of quadrature points of (size = n x n x 2) (optional, adding wb and xb accelerates the integration)

    procedure(myfunc) :: func
    integer, intent(in)  :: d
    character(len=9), dimension(:), allocatable, intent(in) :: t
    integer,  dimension(:), allocatable, intent(in) :: q
    double precision     ,  dimension(:), allocatable, intent(in) ::  p
    double precision     ,  dimension(:), allocatable, intent(in) ::  pmu
    double precision     ,  dimension(:,:), allocatable, intent(in) ::  psigma
    double precision     ,  dimension(:), allocatable, intent(in) ::  r1,r2
    logical , intent(in) :: adap
    integer  ,  dimension(:), allocatable, intent(in) ::  Ni
    double precision, intent(in) :: tol
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: wb    
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: xb    
    double precision     ,  dimension(:), allocatable ::  w,fv,fd
    double precision     ,  dimension(:,:), allocatable ::  x,xx
    integer           :: i,j,Np,maxn
    double precision :: er,s,s0
    integer  ,  dimension(:), allocatable ::  N
    integer  ,  dimension(:,:), allocatable ::  table    
    double precision :: fdmin
    

    er = 100.0d0
    maxn = 100
    allocate(N(d))
    allocate(fd(d))
    do j=1,d
        N(j) = Ni(j)
    enddo

    if (present(xb).and.present(wb)) then
        call quadcreate_aniso_full(q,N,x,w,table,xb,wb)
    else
        call quadcreate_aniso_full(q,N,x,w,table)
    end if

    call quad_scale(q,w,x,r1,r2,t,xx,pmu,psigma)
    Np= size(x,2)

    s = 0.0d0
    do i=1,Np
        s = s + func(xx(:,i),p)*w(i)
    enddo
    
    s0=s
    if (adap) then
        deallocate(xx)
        N=N+2
        do while ((maxval(N,N.gt.1).lt.maxn).and.(er.ge.tol))        

            if (present(xb).and.present(wb)) then
              call quadcreate_aniso_full(q,N,x,w,table,xb,wb)
            else
              call quadcreate_aniso_full(q,N,x,w,table)
            end if

            call quad_scale(q,w,x,r1,r2,t,xx,pmu,psigma)
            Np= size(x,2)

            allocate(fv(Np))
            s = 0.0d0
            do i=1,Np
              fv(i)=func(xx(:,i),p)
              s = s + fv(i)*w(i)
            enddo

            er = dabs((s-s0)/s0)
            do j=1,d
              fd(j) = getmeandiff(fv,x,j,table)
            enddo
            fdmin=minval(fd,fd.ge.0.0d0)
            fd=fd/fdmin
            do j=1,d
              if (fd(j).le.1.05d0) then
                N(j)=N(j)+2
              else if ((fd(j).gt.1.05d0).and.(fd(j).le.1.25d0)) then
                N(j)=N(j)+4
              else if (fd(j).gt.1.25d0) then
                N(j)=N(j)+6
            end if
            enddo

            deallocate(x)
            deallocate(xx)
            deallocate(fv)
            deallocate(w)
            deallocate(table)
            s0 = s
        end do
    end if 
    
    do j=1,d
        if ((t(j)=='lognormal').or.(t(j)=='normal')) then
            s=s/dsqrt((3.1415926535897932d0*2.0d0))
        end if
    enddo
    quad_integral = s    
end function quad_integral

! Function that integrates a multiD function f according to a given distribution
! The function to be integrated has input parameters p
double precision function quad_integral_PDF(func,d,q,p,r1,r2,adap,Ni,tol,wb,xb)
    !implicit none
    ! Using the interface block:
    procedure(myfunc) :: func
    double precision     ,  dimension(:), allocatable ::  w,fv,fd
    double precision     ,  dimension(:), allocatable, intent(in) ::  p
    double precision     ,  dimension(:,:), allocatable ::  x
    integer  ,  dimension(:,:), allocatable ::  table
    integer           :: i,maxn,Np,j
    integer,  dimension(:), allocatable, intent(in) :: q
    double precision     ,  dimension(:), allocatable, intent(in) ::  r1,r2
    logical , intent(in) :: adap
    integer  ,  dimension(:), allocatable, intent(in) ::  Ni
    double precision, intent(in) :: tol
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: wb    
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: xb 
    integer, intent(in)  :: d
    double precision :: er,s,s0
    double precision :: fdmin
    integer  ,  dimension(:), allocatable ::  N

    er = 100.0d0
    maxn = 100
    allocate(N(d))
    allocate(fd(d))
    do j=1,d
        N(j) = Ni(j)
    enddo

    if (present(xb).and.present(wb)) then
        call quadcreate_aniso_full(q,N,x,w,table,xb,wb)
    else
        call quadcreate_aniso_full(q,N,x,w,table)
    end if 
    call quad_scale(q,w,x,r1,r2)

    Np= size(w)
    s = 0.0d0
    do i=1,Np
        s = s + func(x(:,i),p)*w(i)
    enddo
        
    s0=s    
    if (adap) then
        N=N+2
        do while ((maxval(N,N.gt.1).lt.maxn).and.(er.ge.tol))
            if (present(xb).and.present(wb)) then
                call quadcreate_aniso_full(q,N,x,w,table,xb,wb)
            else
                call quadcreate_aniso_full(q,N,x,w,table)
            end if 
            call quad_scale(q,w,x,r1,r2)
            Np= size(w)   
            allocate(fv(Np))

            s = 0.0d0
            do i=1,Np
                fv(i)=func(x(:,i),p)
                s = s + fv(i)*w(i)
            enddo

            er = abs((s-s0)/s0)
            do j=1,d
                fd(j) = getmeandiff(fv,x,j,table)
            enddo
            fdmin=minval(fd,fd.ge.0.0d0)
            fd=fd/fdmin
            do j=1,d
                if (fd(j).le.1.05d0) then
                    N(j)=N(j)+2
                else if ((fd(j).gt.1.05d0).and.(fd(j).le.1.25d0)) then
                    N(j)=N(j)+4
                else if (fd(j).gt.1.25d0) then
                    N(j)=N(j)+6
                end if
            enddo
            deallocate(x)
            deallocate(fv)
            deallocate(w)
            deallocate(table)
            s0 = s
        end do
    end if
    quad_integral_PDF = s
end function quad_integral_PDF

! Function that integrates a multiD function f according to a given distribution
! The function to be integrated has input parameters p
function quad_integral_PDFv(func,d,q,p,r1,r2,adap,Ni,tol,wb,xb)
    !implicit none
    ! Using the interface block:
    procedure(myfuncd) :: func
    double precision     ,  dimension(:), allocatable ::  w
    double precision     ,  dimension(:), allocatable, intent(in) ::  p
    double precision     ,  dimension(:,:), allocatable ::  x
    integer  ,  dimension(:,:), allocatable ::  table
    integer           :: i,maxn,Np,j,dout
    integer,  dimension(:), allocatable, intent(in) :: q
    double precision     ,  dimension(:), allocatable, intent(in) ::  r1,r2
    logical , intent(in) :: adap
    integer  ,  dimension(:), allocatable, intent(in) ::  Ni
    double precision, intent(in) :: tol
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: wb    
    double precision  ,  dimension(:,:,:), allocatable, intent(in), optional :: xb 
    integer, intent(in)  :: d
    double precision :: er
    double precision :: fdmin
    integer  ,  dimension(:), allocatable ::  N
    double precision     ,  dimension(:), allocatable :: s,s0,erv,fd,so
    double precision, dimension(:), allocatable :: quad_integral_PDFv
    double precision     ,  dimension(:,:), allocatable :: fdm,fv

    er = 100.0d0
    maxn = 100
    allocate(N(d))
    do j=1,d
        N(j) = Ni(j)
    enddo

    if (present(xb).and.present(wb)) then
        call quadcreate_aniso_full(q,N,x,w,table,xb,wb)
    else
        call quadcreate_aniso_full(q,N,x,w,table)
    end if 
    call quad_scale(q,w,x,r1,r2)

    Np= size(w)
    s = func(x(:,1),p)*w(1)

    do i=2,Np
        so = func(x(:,i),p)*w(i)
        s = s + so
        deallocate(so)
    enddo    
    
    dout=size(s)
    allocate(s0(dout))
    allocate(erv(dout))
    allocate(fdm(dout,d))
    allocate(fd(d))

    s0=s    
    if (adap) then
        N=N+2
        do while ((maxval(N,N.gt.1).lt.maxn).and.(er.ge.tol))
            if (present(xb).and.present(wb)) then
                call quadcreate_aniso_full(q,N,x,w,table,xb,wb)
            else
                call quadcreate_aniso_full(q,N,x,w,table)
            end if 
            call quad_scale(q,w,x,r1,r2)
            Np= size(w)   
            allocate(fv(dout,Np))

            s = 0.0d0
            do i=1,Np
                  so = func(x(:,i),p)
                  fv(:,i) = so
                  s = s + fv(:,i)*w(i)
                  deallocate(so)
            enddo

            do i=1,dout
                erv(i) = dabs((s(i)-s0(i))/s0(i))
            enddo
            er=maxval(erv,erv.ge.0.0d0)

            do i=1,dout
                do j=1,d
                    fdm(i,j) = getmeandiff(fv(i,:),x,j,table)
                enddo
            enddo

            fd=0.0d0
            do i=1,dout
                fdmin=minval(fdm(i,:),fdm(i,:).ge.0.0d0)
                fdm(i,:)=fdm(i,:)/fdmin
                fd=fd+fdm(i,:)
            enddo
            fd=fd/dble(dout)

            do j=1,d
                if (fd(j).le.1.05d0) then
                    N(j)=N(j)+2
                else if ((fd(j).gt.1.05d0).and.(fd(j).le.1.25d0)) then
                    N(j)=N(j)+4
                else if (fd(j).gt.1.25d0) then
                        N(j)=N(j)+6
                end if
            enddo
            deallocate(x)
            deallocate(fv)
            deallocate(w)
            deallocate(table)
            s0 = s
        end do
    end if

    allocate(quad_integral_PDFv(dout))
    quad_integral_PDFv = s
    
  end function quad_integral_PDFv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! UTILITY FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Utility function for quad_integral_PDF
subroutine quad_scale(tq,w,x,r1,r2,t,xx,pmu,psigma)
    integer, dimension(:), allocatable, intent(in) :: tq
    double precision     ,  dimension(:), allocatable, intent(in) ::  r1,r2
    double precision     ,  dimension(:,:), allocatable, intent(inout) ::  x
    double precision     ,  dimension(:), allocatable, intent(inout) ::  w
    double precision     ,  dimension(:), allocatable, intent(in), optional ::  pmu
    double precision     ,  dimension(:,:), allocatable, intent(in), optional ::  psigma
    character(len=9), dimension(:), allocatable, intent(in), optional :: t
    double precision     ,  dimension(:,:), allocatable, intent(out), optional ::  xx
    integer :: Np,d,j,i

    d=size(x,1)
    Np=size(x,2)
    do j=1,d
        if (tq(j)==1) then          ! Legendre
            x(j,:)=(x(j,:)+1.0d0)/2.0d0*(r2(j)-r1(j))+r1(j)
            w=w*(r2(j)-r1(j))
        else if (tq(j)==2) then     ! Hermite
            x(j,:)=x(j,:)*sqrt(2.0d0)
            w=w*sqrt(2.0d0)
        end if
    enddo

    if (present(t).and.present(xx).and.present(pmu).and.present(psigma)) then
        allocate(xx(d,Np))
        xx=matmul(psigma,x)

        do i=1,Np
            do j=1,d
                if (t(j)=='lognormal') then
                    xx(j,i)=dexp(xx(j,i)+pmu(j))
                    w(i)=w(i)*dexp(-x(j,i)**2/2.0d0)
                else if (t(j)=='normal') then
                    xx(j,i)=xx(j,i)+pmu(j)
                    w(i)=w(i)*dexp(-x(j,i)**2/2.0d0)
              end if
            enddo
        enddo
    end if
end subroutine quad_scale

! Function that returns the mean gradient of a function in a given
! direction
double precision function getmeandiff(fv,x,d,table) 
    double precision     ,  dimension(:,:), allocatable, intent(in) ::  x
    double precision     ,  dimension(:),  intent(in) ::  fv
    integer  ,  dimension(:,:), allocatable, intent(in) ::  table
    integer, intent(in)  :: d
    integer           :: i,Np,j,n,Nd
    double precision  :: s

    Np= size(x,2)
    Nd= size(x,1)
    if (d>Nd) then
        stop
    endif
    s=0.0d0
    n=1
    do j=2,d
        n=n*table(j-1,Np)
    enddo

    do i=1,Np
        if (table(d,i)/=table(d,Np)) then
            s = s + dabs((fv(i+n)-fv(i))/(x(d,i+n)-x(d,i)))
        end if
    enddo

    getmeandiff = s/dble(Np-Np/table(d,Np))
  
end function getmeandiff

end module

