
program cavity_D2Q9
    implicit none
    integer , parameter :: n=50,m=50
    doubleprecision :: f(0:8,0:n,0:m)
    doubleprecision :: feq(0:8,0:n,0:m),rho(0:n,0:m),rho_old(0:n,0:m)
    doubleprecision :: u(0:n,0:m), v(0:n,0:m),u_old(0:n,0:m),v_old(0:n,0:m)
    real :: X(0:n,0:m) , Y(0:n,0:m)
    real :: w(0:8), cx(0:8),cy(0:8)
    real :: alpha , dx , dy , dt , uo , Re , omega , rhoo , eps = 1e-7 &
            & , err_rho , err_u , err_v , PI , Cs
    integer :: i ,  j , kk , mstep


    PI = 4*atan(1.d0)

    err_rho = 0.0
    err_u   = 0.0
    err_v   = 0.0

    !   6 2 5
    !    \|/
    !   3-0-1
    !    /|\
    !   7 4 8

    uo=0.1
    rhoo=5.00
    dx=1.0
    dy=dx
    dt=1.0
    Re=100.0
    Cs = 1/sqrt(3.0)
    alpha=uo*m/Re


    print *, "Re = ", Re , "alpha : " , alpha
    omega=1.0/(3.*alpha+0.5)
    mstep=30000

    w(0)=4./9.
    do i=1,4
        w(i)=1./9.
    end do
    do i=5,8
        w(i)=1./36.
    end do

    ! compute the direction of lattice velocities

    cx(0)=0
    cx(1)=1
    cx(2)=0
    cx(3)=-1
    cx(4)=0
    cx(5)=1
    cx(6)=-1
    cx(7)=-1
    cx(8)=1
    cy(0)=0
    cy(1)=0
    cy(2)=1
    cy(3)=0
    cy(4)=-1
    cy(5)=1
    cy(6)=1
    cy(7)=-1
    cy(8)=-1


    ! initialize rho & u & v
    do j=0,m
        do i=0,n
            rho(i,j)=rhoo
            u(i,j)=0.0
            v(i,j)=0.0
        end do
    end do


    do i=1,n-1
        u(i,m)=uo
        v(i,m)=0.0
    end do

    CALL mesh(n,m,X,Y,dx,dy)

! main loop

    do kk=1,mstep

        rho_old = rho
        u_old   = u
        v_old   = v

        call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
        call streaming(f,n,m)
        ! ——————————–
        call boundary(f,n,m,uo)
        call macro(f,rho,u,v,cx,cy,n,m)

        ! convergence checking

        err_rho = maxval(abs(rho -rho_old))
        err_u   = maxval(abs(u - u_old))
        err_v   = maxval(abs(v - v_old))


        if (mod(kk,1000)==0) then
            print * , "========================================================="
            print *, "the rho residual : " , err_rho
            print *, "the u residual : " , err_u
            print *, "the v residual : " , err_v
            print * , "========================================================="
        end if

        if (err_rho .LE. eps .AND. err_u .LE. eps .AND. err_v .LE. eps) then
            print * , "the solution has been converged at iteration : " , kk
            exit
        end if
    END DO
! end of the main loop

    CALL output(n,m,dx,dy,rho,u,v,X,Y)
    CALL section(n,m,u,X,Y,uo)

    contains

! end of the main program
    subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
        integer , INTENT(IN) :: n , m
        doubleprecision , INTENT(INOUT) :: f(0:8,0:n,0:m) , feq(0:8,0:n,0:m)
        doubleprecision , INTENT(INOUT) :: rho(0:n,0:m) , u(0:n,0:m), v(0:n,0:m)
        real , INTENT(IN) :: w(0:8), cx(0:8),cy(0:8)
        real , INTENT(IN) :: omega
        integer :: i, j , k
        real :: t1 , t2

            DO i=0,n
                DO j=0,m
                    t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
                    DO k=0,8
                        t2=u(i,j)*cx(k)+v(i,j)*cy(k)
                        feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
                        f(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
                    END DO
                END DO
            END DO
    end subroutine collesion


    subroutine streaming(f,n,m)
        integer , INTENT(IN) :: n , m
        doubleprecision , INTENT(INOUT) :: f(0:8,0:n,0:m)
        integer :: i , j

        ! streaming
        DO j=0,m
            DO i=n,1,-1 !RIGHT TO LEFT
                f(1,i,j)=f(1,i-1,j)
            END DO
            DO i=0,n-1 !LEFT TO RIGHT
                f(3,i,j)=f(3,i+1,j)
            END DO
        END DO

        DO j=m,1,-1 !TOP TO BOTTOM
            DO i=0,n
                f(2,i,j)=f(2,i,j-1)
            END DO
            DO i=n,1,-1
                f(5,i,j)=f(5,i-1,j-1)
            END DO
            DO i=0,n-1
                f(6,i,j)=f(6,i+1,j-1)
            END DO
        END DO

        DO j=0,m-1 !BOTTOM TO TOP
            DO i=0,n
                f(4,i,j)=f(4,i,j+1)
            END DO
            DO i=0,n-1
                f(7,i,j)=f(7,i+1,j+1)
            END DO
            DO i=n,1,-1
                f(8,i,j)=f(8,i-1,j+1)
            END DO
        END DO

    end subroutine streaming

    subroutine boundary(f,n,m,uo)
        doubleprecision , INTENT(INOUT) :: f(0:8,0:n,0:m)
        integer , INTENT(IN) :: n , m
        real , INTENT (IN)   :: uo
        integer :: i, j , k
        real :: rhon , t1 , t2

        do j=0,m
        ! bounce back on west boundary
            f(1,0,j)=f(3,0,j)
            f(5,0,j)=f(7,0,j)
            f(8,0,j)=f(6,0,j)
        ! bounce back on east boundary
            f(3,n,j)=f(1,n,j)
            f(7,n,j)=f(5,n,j)
            f(6,n,j)=f(8,n,j)
        end do
        ! bounce back on south boundary
        do i=0,n
            f(2,i,0)=f(4,i,0)
            f(5,i,0)=f(7,i,0)
            f(6,i,0)=f(8,i,0)
        end do
        ! moving lid, north boundary
        do i=1,n-1
            rhon=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
            f(4,i,m)=f(2,i,m)
            f(8,i,m)=f(6,i,m)+rhon*uo/6.0
            f(7,i,m)=f(5,i,m)-rhon*uo/6.0
        end do

    end subroutine boundary


    subroutine macro(f,rho,u,v,cx,cy,n,m)
        integer , INTENT(IN) :: m , n
        doubleprecision , INTENT(IN) :: f(0:8,0:n,0:m)
        doubleprecision , INTENT(INOUT) :: rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m)
        real , INTENT(IN) :: cx(0:8),cy(0:8)
        integer :: i , j , k
        real :: ssum , usum , vsum

        do j=0,m
            do i=0,n
                ssum=0.0
            do k=0,8
                ssum=ssum+f(k,i,j)
            end do
                rho(i,j)=ssum
            end do
        end do

        DO i=1,n
            DO j=1,m-1
                usum=0.0
                vsum=0.0
                DO k=0,8
                    usum=usum+f(k,i,j)*cx(k)
                    vsum=vsum+f(k,i,j)*cy(k)
                END DO
                u(i,j)=usum/rho(i,j)
                v(i,j)=vsum/rho(i,j)
            END DO
        END DO

    end subroutine macro

    subroutine mesh(n,m,X,Y,dx,dy)
        integer , INTENT(IN) :: m , n
        real , INTENT(IN) :: dx , dy
        real ,INTENT(OUT) :: X(0:n,0:m) , Y(0:n,0:m)
        integer :: i,j

        do i=0,n
            do j=0,m
                X(i,j) = real(i)*dx
                Y(i,j) = real(j)*dy
            end do
        end do

    end subroutine mesh

    subroutine output(n,m,dx,dy,rho,u,v,X,Y)
        integer , INTENT(IN) :: m , n
        real , INTENT(IN) :: dx , dy
        doubleprecision , INTENT(IN) :: rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m)
        real , INTENT(IN):: X(0:n,0:m) , Y(0:n,0:m)
        integer :: i,j

        open(unit=1,file="rho.dat")

        write(1,*) "VARIABLE = X , Y , U , V"
        write(1,*) "ZONE " , ", I=" , n+1 , ", J=" , m+1

        do i=0,n
            do j=0,m
                write(1,*) X(i,j) , Y(i,j) , u(i,j) , v(i,j)
            end do
        end do


    end subroutine output

    subroutine section(n,m,u,X,Y,uo)
        integer,INTENT(IN) :: n,m
        doubleprecision,INTENT(IN)::u(0:n,0:m)
        real,INTENT(IN):: X(0:n,0:m) , Y(0:n,0:m),uo

        open(unit=1,file="varticalSEC.dat") ! x=0.5
        write(1,*) "VARIABLE Y , U"

        do j=0,m
            write(1,*) Y(n/2,j)/real(m) , u(n/2,j)/uo
        end do

        open(unit=2,file="horizontalSEC.dat") ! y=0.5
        write(2,*) "VARIABLE X , U"

        do i=0,n
            write(2,*) X(i,m/2)/real(n) , u(i,m/2)/uo
        end do

    end subroutine

end program cavity_D2Q9

