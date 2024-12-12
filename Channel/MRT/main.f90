
program channel_MRT
    implicit none
    integer , parameter :: n=2000,m=80
    doubleprecision :: f(0:8,0:n,0:m),feq(0:8,0:n,0:m)
    doubleprecision :: rho(0:n,0:m),rho_old(0:n,0:m)
    doubleprecision :: u(0:n,0:m),v(0:n,0:m),u_old(0:n,0:m),v_old(0:n,0:m),L2
    real :: X(0:n,0:m),Y(0:n,0:m),tminv(0:8,0:8),sm(0:8),tm(0:8,0:8),stmiv(0:8,0:8)
    real :: w(0:8),cx(0:8),cy(0:8),ev(0:8,0:8),start,finish
    real :: nu , dx , dy , dt , uo , Re , omega , rhoo , eps = 1e-8 &
            & , err_rho , err_u , err_v , a1 , sumcc,tau,PI ,Cs
    integer :: i,j,k, kk , mstep

    PI = 4*atan(1.d0)


            tm(0,:) = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
            tm(1,:) = (/-4.0,-1.0,-1.0,-1.0,-1.0,2.0,2.0,2.0,2.0/)
            tm(2,:) = (/4.0,-2.0,-2.0,-2.0,-2.0,1.0,1.0,1.0,1.0/)
            tm(3,:) = (/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
            tm(4,:) = (/0.0,-2.0,0.0,2.0,0.0,1.0,-1.0,-1.0,1.0/)
            tm(5,:) = (/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
            tm(6,:) = (/0.0,0.0,-2.0,0.0,2.0,1.0,1.0,-1.0,-1.0/)
            tm(7,:) = (/0.0,1.0,-1.0,1.0,-1.0,0.0,0.0,0.0,0.0/)
            tm(8,:) = (/0.0,0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0/)

            a1 = 1.0/36.0

            tminv(0,:) = (/4.0*a1,-4.0*a1,4.0*a1,0.0,0.0,0.0,0.0,0.0,0.0/)
            tminv(1,:) = (/4.0*a1,-a1,-2.0*a1,6.0*a1,-6.0*a1,0.0,0.0,9.0*a1,0.0/)
            tminv(2,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,6.0*a1,-6.0*a1,-9.0*a1,0.0/)
            tminv(3,:) = (/4.0*a1,-a1,-2.0*a1,-6.0*a1,6.0*a1,0.0,0.0,9.0*a1,0.0/)
            tminv(4,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,-6.0*a1,6.0*a1,-9.0*a1,0.0/)
            tminv(5,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,6.0*a1,3.0*a1,0.0,9.0*a1/)
            tminv(6,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,6.0*a1,3.0*a1,0.0,-9.0*a1/)
            tminv(7,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,-6.0*a1,-3.0*a1,0.0,9.0*a1/)
            tminv(8,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,-6.0*a1,-3.0*a1,0.0,-9.0*a1/)


    err_rho = 0.0
    err_u   = 0.0
    err_v   = 0.0
    mstep=50000

    !   6 2 5
    !    \|/
    !   3-0-1
    !    /|\
    !   7 4 8

    w(0)=4./9.
    do i=1,4
        w(i)=1./9.
    end do
    do i=5,8
        w(i)=1./36.
    end do

    ! compute the direction of lattice velocities
    cx(:) = (/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
    cy(:) = (/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)


    CALL init(n,m,u,v,rho,omega,tau,dx,dy,uo,nu)
    CALL mesh(n,m,X,Y,dx,dy)
    CALL MRT(n,m,tminv,sm,tm,stmiv,tau,nu)

    open(unit=10,file='error.dat')
    write(10,*) "VARIABLE =  step , error"

! main loop
    do kk=1,mstep

        rho_old = rho
        u_old   = u
        v_old   = v

        call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,tm,tminv,stmiv)
        call streaming(f,n,m)
        ! ——————————–
        call boundary(f,n,m,uo)
        call macro(f,rho,u,v,cx,cy,n,m)

        ! convergence checking

        err_rho = maxval(abs(rho - rho_old))
        err_u   = maxval(abs(u   - u_old))
        err_v   = maxval(abs(v   - v_old))

        CALL CPU_TIME(start)

        if (mod(kk,100)==0) then
            print * , "========================================================="
            print *, "the rho residual : " , err_rho
            print *, "the u residual   : " , err_u
            print *, "the v residual   : " , err_v
            print *, "time step        : " , kk
            CALL check_convergence(n,m,kk,u,u_old,v,v_old,L2)
            if (L2 .LT. eps) then
                print * , char(7)
                print * , "the solution has been converged at iteration : " , kk
                exit
            end if
            print *, "the relative error :  " , L2
            print * , "========================================================="
        end if

        if(kk==mstep) print * , char(7)

    END DO
    CALL CPU_TIME(finish)
! end of the main loop

    CALL output(n,m,dx,dy,rho,u,v,X,Y)

    CALL analytical(n,m,X,Y,uo)

    print * , "the cpu time : " , finish-start

    contains

! end of the main program
    subroutine init(n,m,u,v,rho,omega,tau,dx,dy,uo,nu)
        integer,INTENT(IN)::n,m
        doubleprecision,INTENT(INOUT)::u(0:n,0:m),v(0:n,0:m),rho(0:n,0:m)
        real,INTENT(INOUT):: omega,tau,dx,dy,uo,nu
        real:: rhoo,dt,Re
        integer::i,j

        uo=0.1
        rhoo=1.00
        dx=1.0
        dy=dx
        dt=1.0
        Re = 100.0
        nu=uo*m/Re


        print *, "Re = ", Re , "nu = " , nu
        omega=1.0/(3.*nu+0.5)
        tau  = 1.0/omega

        do j=0,m
            do i=0,n
                rho(i,j)=rhoo
                u(i,j)  =0.0
                v(i,j)  =0.0
            end do
        end do

        do j=1,m-1
            u(0,j)=uo
            v(0,j)=0.0
        end do

    end subroutine

    subroutine MRT(n,m,tminv,sm,tm,stmiv,tau,nu)
        integer,INTENT(IN) ::n,m
        real,INTENT(INOUT) :: tminv(0:8,0:8),sm(0:8),tm(0:8,0:8),stmiv(0:8,0:8),tau,nu
        real::s3,s5,s7,s8,dt=1.0
        integer:: i,j,k,l

            Cs = 1.0/sqrt(3.0)

            s3 = 1.0
            s5 = s3
            s7 = 1.0/tau
            s8 = s7

            sm(:) = (/0.8,1.,1.,s3,1.2,s5,1.2,s7,s8/)

            do i=0,8
                do j=0,8
                    stmiv(i,j) = tminv(i,j)*sm(j)
                end do
            end do


            do i=0,8
                do j=0,8
                    sumcc = 0.0
                    do k=0,8
                        sumcc = sumcc+tminv(i,k)*tm(k,j)
                    end do
                    ev(i,j) = sumcc
                end do
            end do

    end subroutine

    subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,tm,tminv,stmiv)
            integer , INTENT(IN) :: n , m
            doubleprecision , INTENT(INOUT) :: f(0:8,0:n,0:m),feq(0:8,0:n,0:m),rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m)
            doubleprecision :: fmom(0:8,0:n,0:m),fmeq(0:8,0:n,0:m)
            real , INTENT(IN) :: w(0:8),cx(0:8),cy(0:8),omega,tminv(0:8,0:8),tm(0:8,0:8),stmiv(0:8,0:8)
            real :: suma,sumb
            integer :: i,j,k,l

            ! calculate the feq moments:
            do i=0,n
                do j=0,m
                    fmeq(0,i,j) = rho(i,j)
                    fmeq(1,i,j) = rho(i,j)*(-2.0+3.0*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
                    fmeq(2,i,j) = rho(i,j)*(1.0-3.0*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
                    fmeq(3,i,j) = rho(i,j)*u(i,j)
                    fmeq(4,i,j) = -rho(i,j)*u(i,j)
                    fmeq(5,i,j) = rho(i,j)*v(i,j)
                    fmeq(6,i,j) = -rho(i,j)*v(i,j)
                    fmeq(7,i,j) = rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
                    fmeq(8,i,j) = rho(i,j)*u(i,j)*v(i,j)
                end do
            end do

            ! calculate the moment:
            do i=0,n
                do j=0,m
                    do k=0,8
                        suma  = 0.0
                        do l=0,8
                            suma = suma + tm(k,l)*f(l,i,j)
                        end do
                        fmom(k,i,j) = suma
                    end do
                end do
            end do

            ! collision

            do i=0,n
                do j=0,m
                    do k=0,8
                        sumb = 0.0
                        do l=0,8
                            sumb = sumb +stmiv(k,l)*(fmom(l,i,j)-fmeq(l,i,j))
                        end do
                        f(k,i,j) = f(k,i,j)-sumb
                    end do
                end do
            end do

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
        real :: rhow , t1 , t2

        ! bounce back on west
        do j=0,m
            rhow = (f(0,0,j)+f(2,0,j)+f(4,0,j)+2.0*(f(3,0,j)+f(6,0,j)+f(7,0,j)))/(1.0-uo)
            f(1,0,j) = f(3,0,j)+2.0*rhow*uo/3.0
            f(5,0,j) = f(7,0,j)-0.5*(f(2,i,j)-f(4,i,j))+rhow*uo/6.0
            f(8,0,j) = f(6,0,j)+0.5*(f(2,i,j)-f(4,i,j))+rhow*uo/6.0
        end do

        ! bounce back on south boundary
        do i=0,n
            f(2,i,0)=f(4,i,0)
            f(5,i,0)=f(7,i,0)
            f(6,i,0)=f(8,i,0)
        end do

        ! bounce back on north boundary
        do i=0,n
            f(4,i,m) = f(2,i,m)
            f(8,i,m) = f(6,i,m)
            f(7,i,m) = f(5,i,m)
        end do

        ! extrapolation on outer region boundary

        do j=1,m
            f(1,n,j) = 2.0*f(1,n-1,j)-f(1,n-2,j)
            f(5,n,j) = 2.0*f(5,n-1,j)-f(5,n-2,j)
            f(8,n,j) = 2.0*f(8,n-1,j)-f(8,n-2,j)
        end do

        do j=0,m
            v(n,j) = 0.0
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


        do j=1,m
            v(n,j) = 0.0
        end do

        do i=0,n
            u(i,0) = 0.0
            u(i,m) = 0.0
        end do

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

    subroutine check_convergence(n,m,kk,u1,u2,v1,v2,L2)
        implicit none
        integer , INTENT(IN) :: n,m,kk
        doubleprecision , INTENT(IN) :: u1(0:n,0:m),u2(0:n,0:m),v1(0:n,0:m),v2(0:n,0:m)
        doubleprecision , INTENT(OUT) :: L2
        real :: tmp1,tmp2
        integer :: i,j

        tmp1 = 0.0
        tmp2 = 0.0

        do i=0,n
            do j=0,m
                tmp1 = tmp1 + (u1(i,j)-u2(i,j))**2+(v1(i,j)-v2(i,j))**2
                tmp2 = tmp2 + u2(i,j)**2+v2(i,j)**2
            end do
        end do

        if(tmp2==0.0) then
            L2 = 1.0
        else
            L2 = sqrt(tmp1/tmp2)
        endif

        write(10,*) kk , L2

    end subroutine


    subroutine output(n,m,dx,dy,rho,u,v,X,Y)
        integer , INTENT(IN) :: m , n
        real , INTENT(IN) :: dx , dy
        doubleprecision , INTENT(IN) :: rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m)
        real , INTENT(IN):: X(0:n,0:m) , Y(0:n,0:m)
        integer :: i,j

        open(unit=1,file="output.dat")

        write(1, '(A)') 'TITLE = "output"'
        write(1, '(A)') 'VARIABLES = "X", "Y" , "U" , "V"'

        write(1, '(A, I0, A, I0, A)') 'ZONE I=', n+1, ', J=', m+1, ', F=POINT'

        do j=0,m
            do i=0,n
                write(1,*) X(i,j) , Y(i,j) , u(i,j) , v(i,j)
            end do
        end do

    end subroutine output

    subroutine analytical(n,m,X,Y,uo)
        integer,INTENT(IN)::n,m
        real,INTENT(IN):: X(0:n,0:m),Y(0:n,0:m),uo
        real:: H,Uin,yi

        Uin = 0.02

        H  = 1.0/2.0

        open(unit=1,file="section.dat") ! at x = 0.5
        write(1,*) "VARIABLE Y , U"
        write(1,*) "ZONE"

        do j=0,m
            write(1,*) (Y(int(n/6),j)/real(m)) , u(int(n/6),j)/maxval(u(int(n/6),:))
        end do

        write(1,*) "ZONE"

        do j=0,m
            write(1,*) (Y(int(2*n/6),j)/real(m)) , u(int(2*n/6),j)/maxval(u(int(2*n/6),:))
        end do

        write(1,*) "ZONE"

        do j=0,m
            write(1,*) (Y(int(3*n/6),j)/real(m)) , u(int(3*n/6),j)/maxval(u(int(3*n/6),:))
        end do

        write(1,*) "ZONE"

        do j=0,m
            write(1,*) (Y(int(4*n/6),j)/real(m)) , u(int(4*n/6),j)/maxval(u(int(4*n/6),:))
        end do

        write(1,*) "ZONE"

        do j=0,m
            write(1,*) (Y(int(5*n/6),j)/real(m)) , u(int(5*n/6),j)/maxval(u(int(5*n/6),:))
        end do

        write(1,*) "ZONE"

        do j=0,m
            yi = (Y(n/2,j)/real(m)-0.5)
            write(1,*) Y(n/2,j)/real(m) , (1-(yi**2)/H**2)
        end do

    end subroutine

end program channel_MRT

